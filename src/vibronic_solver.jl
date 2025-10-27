#
#~~~~~~~~~~~~~~~~~~~~~~~~~~ VIBRONIC SOLVER FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~#
#                       SINC-DVR FOR CONTRACTED SOLUTION                       #
#
## sin(x)/x for sinc-DVR
function sinc_function(x::Float64, xj::Float64, h::Float64)::Float64
    # Avoid division by zero: handle x == x_j separately
    if abs(x - xj) < 1e-10
        return 1.0
    else
        return sin(pi * (x - xj) / h) / (pi * (x - xj) / h)
    end
end
#
## compute the uncoupled (no NAC) kinetic energy matrix in the sinc-DVR basis: T
function SincDVR_KinMat(vmax::Int,r_min::Float64, r_max::Float64, factor::Float64)::Matrix{Float64}
    #
    T = zeros(Float64, vmax, vmax)
    #
    ## compute quadrature points
    R = LinRange(r_min, r_max,vmax)
    #
    h = R[2] - R[1]
    #
    @inbounds for i=1:vmax
        for j=i:vmax
            if i == j
                T[i,i] = pi^2 / (3 * h^2)
            else
                T[i,j] = 2 * ( (-1)^(i-j) / (i-j)^2 ) / h^2
                T[j,i] = 2 * ( (-1)^(i-j) / (i-j)^2 ) / h^2
            end
        end
    end
    #
    return factor * T
end
#
## compute the potential energy matrix in the sinc-DVR basis: V
function SincDVR_PotMat(vmax::Int, fPot, r_min::Float64, r_max::Float64)::Diagonal(Float64,Vector{Float64})
    #
    ## evaluate quadrature points
    R = LinRange(r_min,r_max,vmax)
    #
    ## compute potential on the quadrature points
    V = Diagonal(fPot.(R))
    #
    return V
end
#
## diagonalise the total sinc-DVR Hamiltonian: H = T + V
function solve_vibrational_schrodinger!( T                :: Matrix{Float64}, 
                                         V                :: Diagonal{Float64, Vector{Float64}}, 
                                         r                :: Vector{Float64}, 
                                         h                :: Float64, 
                                         electronic_state :: Int64,
                                         φ                :: Vector{Vector{Float64}},
                                         state_v_wfn      :: Array{Float64, 3},
                                         state_v_ener     :: Matrix{Float64}
                                         )                #:: Tuple{Vector{Float64}, Matrix{Float64}}
    """
    Solve the vibrational Schrödinger equation Hψ = Eψ using DVR.
    
    Parameters:
        T: Kinetic energy matrix (square matrix)
        V: Potential energy matrix (diagonal matrix, same size as T)
    
    Returns:
        E: Eigenvalues (vibrational energy levels)
        Ψ: Eigenvectors (vibrational wavefunctions in DVR basis)
    """
    #
    ## Construct the total Hamiltonian
    H = T + V
    #
    ## Diagonalize the Hamiltonian
    eigenvalues, eigenvectors = eigen(H)
    #
    ## number of eigenvalues
    Nv = length(eigenvalues)
    #
    ## compute uncoupled vibronic wavefunctions as basis and populate an array
    @inbounds for idx=1:Nv
        #
        v = idx - 1
        #
        state_v_wfn[electronic_state,v+1,:] .= Vibronic_Wavefunction(v, eigenvectors, φ, r)
        #
        state_v_ener[electronic_state,v+1] = eigenvalues[idx]
    end
    #
    # return eigenvalues, eigenvectors
end
#
## normalise wavefunction
function normalise_wavefunction(wfn::Vector{Float64}, r::Vector{Float64})
    norm = sqrt(simps(abs2.(wfn), r[1], r[end]))
    return wfn ./ norm
end
#
## compute the r-dependent vibronic wavefunction from the sinc-DVR basis
function Vibronic_Wavefunction( v                 :: Int64, 
                                eigenvectors      :: Matrix{Float64}, 
                                φ                 :: Vector{Vector{Float64}}, 
                                r                 :: Vector{Float64}; 
                                thresh = 1e-100
                              )                   :: Vector{Float64}
    #
    ## sum the vibronic basis with the weights coeff, only compute if coefficient is above threshold
    Ψv = zeros(Float64, length(r))
    #
    coeff = eigenvectors[:,v+1]
    #
    @inbounds for i in eachindex(coeff)
        Ψv .+= coeff[i] .* φ[i]
    end
    #
    return normalise_wavefunction(Ψv, r)
end
#
## method function in solving the pure vibronic problem for a set of electronic states
function vibronic_eigensolver(r::Vector{Float64}, V::Array{Float64, 3}, vmax::Int64, states::Vector{Int64}, m1::Float64, m2::Float64)
    #
    ## compute number of electronic states
    N_el = length(states)
    #
    ## initialise vibronic wavefunction object
    state_v_wfn = Array{Float64}(undef, N_el, vmax, vmax)
    #
    ## initialise vibronic energy object
    state_v_ener = Array{Float64}(undef, N_el, vmax)
    #
    ## compute the grid seperation
    h = r[2] - r[1]
    #
    ## compute the kinetic energy matrix in sinc-DVR basis
    T = SincDVR_KinMat(vmax, r[1], r[end], KE_factor(m1,m2))
    #
    ## precompute Sinc-DVR basis functions
    φ = [sinc_function.(r, x, h) for x in r]
    #
    ## compute the potential matrix in sinc-DVR basis (V) for all electronic states
    ## and then diagonalise the total Hamiltonian H = T + V
    for state in states
        #
        V_state = Diagonal(V[:,state,state])
        #
        ## solve vibrational Schrodinger equation for the electronic state
        solve_vibrational_schrodinger!(T, V_state, r, h, state, φ, state_v_wfn, state_v_ener)
    end
    #
    ## return energies and vibronic wavefunctions
    return state_v_wfn, state_v_ener
end
#
## visualise the vibronic solution by plotting wavefunctions, energies, and PECs
function plot_vibronic_solution(r::Vector{Float64}, V::Array{Float64, 3}, states::Vector{Int64}, E, psi, vmax::Array{Int64})
    #
    plt.figure()
    for state in states
        #
        plt.plot(r,V[:,state,state])
        #
        for v_idx=1:vmax[state]
            #
            Ev = E[state,v_idx]
            #
            Psi_v = psi[state,v_idx,:]
            #
            ## compute the scaling for the wavefunction to plot
            ∂E = E[state,v_idx+1] - Ev
            #
            scale = 0.4 * ∂E
            #
            prob_density = (Psi_v).^2
            prob_density_norm = prob_density ./ maximum(prob_density)
            #
            ## compute where the wavefunction drops below 1e-4
            mask = prob_density_norm .> 1e-4
            idx_i = findfirst(mask)
            idx_f = lastindex(prob_density_norm) - findfirst(reverse(mask)) + 1
            #
            plt.plot(r[idx_i:idx_f], scale .* prob_density_norm[idx_i:idx_f] .+ Ev)
        end
    end
end
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~ COUPLED VIBRONIC PROBLEM ~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#
## solve the coupled vibronic schrodginer equation. The equation is a matrix
## operator equation, where the matrix elements are in the electronic basis.
#
## initialise the (coupled) kinetic energy hamiltonian:
#                      -->            <--        -->
## T = -hbar^2/2mu * (d2/dr2 + W^2 - (d/dr W - W d/dr)))
#
# function Build_nonadiabatic_T_V() #:: TYPE DECLARE NEEDED
#     #
#     ## flatten multidimensional electronic+vibrational indexing
#     function serialise_electronic_vib_indices(state, v_idx, contraction_array)
#         return ( state - 1 ) * contraction_array[ state ] + v_idx
#     end
#     #
#     ## initialise contracted vibronic basis size
#     vmax = Calculation["vibronic_solver"].contraction
#     contracted_vibronic_dim = sum(vmax)
#     #
#     ## initialise the kinetic energy and electronic Hamiltonians
#     T = zeros(Float64, contracted_vibronic_dim, contracted_vibronic_dim)
#     V = zeros(Float64, contracted_vibronic_dim, contracted_vibronic_dim)
#     #
#     ## initialise the number of electronic states
#     Nstates = length( Calculation["method"].states )
#     #
#     ## initialise the NAC squared matrix: second DDR
#     W2 = zeros(Float64, length(r), Nstates, Nstates)
#     K = map(idx -> NACMat[idx,:,:] * NACMat[idx,:,:], collect(1:lastindex(r)))
#     for i=1:Nstates
#         for j=i:Nstates
#             W2[:,i,j] .= [K[idx][i,j] for idx=1:lastindex(r)]
#         end
#     end
#     #
#     ## precompute wavefunction derivatives
#     wfn_ddr   = zeros(Float64, lastindex(r), contracted_vibronic_dim) 
#     wfn_d2dr2 = zeros(Float64, lastindex(r), contracted_vibronic_dim) 
#     #
#     for i=1:Nstates
#         for v_idx=1:vmax[i]
#             #
#             state_v = serialise_electronic_vib_indices(i, v_idx, vmax)
#             #
#             wfn_ddr[:, state_v] .= FiniteDifference(r, contr_vib_wfn[i,v_idx,:], 1)
#             #
#             wfn_d2dr2[:, state_v] .= FiniteDifference(r, contr_vib_wfn[i,v_idx,:], 2)
#         end
#     end
#     #
#     ## electronically diagonal terms
#     for i=1:Nstates
#         W2_ii = W2[:,i,i]
#         #
#         V_ii = PotMat[:,i,i]
#         #
#         for v_idx=1:vmax[i]
#             for v_jdx=v_idx:vmax[i]
#                 #
#                 row    = serialise_electronic_vib_indices(i, v_idx, vmax)
#                 column = serialise_electronic_vib_indices(i, v_jdx, vmax)
#                 #
#                 d2dr2 = contr_vib_wfn[i,v_idx,:] .* wfn_d2dr2[:, column]
#                 #
#                 W2_int = contr_vib_wfn[i,v_idx,:] .* W2_ii .* contr_vib_wfn[i,v_jdx,:]
#                 #
#                 integrand = d2dr2  .+ W2_int
#                 #
#                 matel =  simps(integrand, r[1], r[end])
#                 #
#                 T[row, column] = matel
#                 T[column, row] = matel
#                 #
#                 ## now for the potential
#                 V_integrand = contr_vib_wfn[i,v_idx,:] .* V_ii .* contr_vib_wfn[i,v_jdx,:]
#                 V_matel = simps(V_integrand, r[1], r[end])
#                 #
#                 V[row, column] = V_matel
#                 V[column, row] = V_matel
#             end
#         end
#     end
#     #
#     ## electronically off-diagonal terms
#     for i=1:Nstates
#         for j=i+1:Nstates
#             #
#             Wij = NACMat[:,i,j]
#             #
#             for v_idx=1:vmax[i]
#                 for v_jdx=1:vmax[j]
#                     #
#                     row    = serialise_electronic_vib_indices(i, v_idx, vmax)
#                     column = serialise_electronic_vib_indices(j, v_jdx, vmax)
#                     #
#                     W2_ij = contr_vib_wfn[i,v_idx,:] .* W2[:,i,j] .* contr_vib_wfn[j,v_jdx,:]
#                     #
#                     dyi_dr = wfn_ddr[:, row]
#                     dyj_dr = wfn_ddr[:, column]
#                     #
#                     ddr = (dyi_dr .* Wij .* contr_vib_wfn[j,v_jdx,:]) .- (contr_vib_wfn[i,v_idx,:] .* Wij .* dyj_dr)
#                     #
#                     integrand = W2_ij .- ddr
#                     #
#                     matel = simps(integrand, r[1], r[end])
#                     #
#                     T[row,column] = matel
#                     T[column,row] = matel
#                 end
#             end
#         end
#     end
#     #
#     return T, V
# end

@inline function flatten_index(state, local_idx, contraction_array)
    return (state - 1) * contraction_array[state] + local_idx
end
#
## flatten multidimensional electronic+vibrational indexing
@inline function serialise_electronic_vib_indices(state, v_idx, contraction_array)
    return ( state - 1 ) * contraction_array[ state ] + v_idx
end
#
## flatten multidimensional total rovibronic Hamiltonian indexing
@inline function serialise_vibronic_rotational_indices(state, v_idx, rot_idx, contraction_array, rot_dims)
    state_v = ( state - 1 ) * contraction_array[ state ] + v_idx
    #
    return ( state_v - 1 ) * rot_dims[ state ] + rot_idx 
end
#
function build_coupled_vibronic_Hamiltonian(contr_vib_wfn) #:: TYPE DECLARE NEEDED
    #
    ## initialise contracted vibronic basis size
    vmax = Calculation["vibronic_solver"].contraction
    contracted_vibronic_dim = sum(vmax)
    #
    ## initialise the kinetic energy and electronic Hamiltonians
    T = zeros(Float64, contracted_vibronic_dim, contracted_vibronic_dim)
    V = zeros(Float64, contracted_vibronic_dim, contracted_vibronic_dim)
    B = zeros(Float64, contracted_vibronic_dim, contracted_vibronic_dim)
    #
    ## initialise the number of electronic states
    Nstates = length( Calculation["method"].states )
    #
    ## initialise the NAC squared matrix: second DDR
    W2 = zeros(Float64, length(r), Nstates, Nstates)
    @views for idx in 1:lastindex(r)
        A = NACMat[idx, :, :]  # avoids copying
        Wsqrd = A * A                     # full matrix multiplication
        for i in 1:Nstates
            for j in i:Nstates
                W2[idx, i, j] = Wsqrd[i, j]
            end
        end
    end
    #
    ## precompute wavefunction derivatives
    wfn_ddr   = zeros(Float64, lastindex(r), contracted_vibronic_dim) 
    wfn_d2dr2 = zeros(Float64, lastindex(r), contracted_vibronic_dim) 
    #
    for i=1:Nstates
        for v_idx=1:vmax[i]
            #
            state_v = serialise_electronic_vib_indices(i, v_idx, vmax)
            #
            wfn_ddr[:, state_v] .= FiniteDifference(r, contr_vib_wfn[i,v_idx,:], 1)
            #
            wfn_d2dr2[:, state_v] .= FiniteDifference(r, contr_vib_wfn[i,v_idx,:], 2)
        end
    end
    #
    ## electronically diagonal terms
    for i=1:Nstates
        W2_ii = W2[:,i,i]
        #
        V_ii = PotMat[:,i,i]
        #
        @views for v_idx=1:vmax[i]
            for v_jdx=v_idx:vmax[i]
                #
                row    = serialise_electronic_vib_indices(i, v_idx, vmax)
                column = serialise_electronic_vib_indices(i, v_jdx, vmax)
                #
                ## kinetic energy
                d2dr2 = contr_vib_wfn[i,v_idx,:] .* wfn_d2dr2[:, column]
                #
                W2_int = contr_vib_wfn[i,v_idx,:] .* W2_ii .* contr_vib_wfn[i,v_jdx,:]
                #
                integrand = d2dr2  .+ W2_int
                #
                matel =  simps(integrand, r[1], r[end])
                #
                T[row, column] = matel
                T[column, row] = matel
                #
                ## now for the potential
                V_integrand = contr_vib_wfn[i,v_idx,:] .* V_ii .* contr_vib_wfn[i,v_jdx,:]
                V_matel = simps(V_integrand, r[1], r[end])
                #
                V[row, column] = V_matel
                V[column, row] = V_matel
                #
                ## now for the rotational constant
                B_matel = simps(contr_vib_wfn[i,v_idx,:] .* r.^(-2) .* contr_vib_wfn[i,v_jdx,:],  r[1], r[end])
                B[row, column] = B_matel
                B[column, row] = B_matel
            end
        end
    end
    #
    ## electronically off-diagonal terms
    for i=1:Nstates
        for j=i+1:Nstates
            #
            Wij = NACMat[:,i,j]
            #
            for v_idx=1:vmax[i]
                for v_jdx=1:vmax[j]
                    #
                    row    = serialise_electronic_vib_indices(i, v_idx, vmax)
                    column = serialise_electronic_vib_indices(j, v_jdx, vmax)
                    #
                    W2_ij = contr_vib_wfn[i,v_idx,:] .* W2[:,i,j] .* contr_vib_wfn[j,v_jdx,:]
                    #
                    dyi_dr = wfn_ddr[:, row]
                    dyj_dr = wfn_ddr[:, column]
                    #
                    ddr = (dyi_dr .* Wij .* contr_vib_wfn[j,v_jdx,:]) .- (contr_vib_wfn[i,v_idx,:] .* Wij .* dyj_dr)
                    #
                    integrand = W2_ij .- ddr
                    #
                    matel = simps(integrand, r[1], r[end])
                    #
                    T[row,column] = matel
                    T[column,row] = matel
                end
            end
        end
    end
    #
    ##
    kefac = -KE_factor(Calculation["method"].atoms...)
    return T * kefac, V, B
end
#
function check_omega_condition(J::Float64, Omega::Float64, Lambda::Float64, S::Float64)::Bool
    """
    Checks if the given Omega value satisfies the rovibronic selection rule:
    |Omega| <= min(J, |Lambda| + S)

    Arguments:
        J: Total angular momentum quantum number.
        Omega: Projection of total angular momentum along the internuclear axis (Lambda + Sigma).
        Lambda: Projection of electronic orbital angular momentum along the internuclear axis.
        S: Total electron spin angular momentum.

    Returns:
        True if the condition is met, False otherwise.
    """
    
    # Calculate the right-hand side of the inequality
    # abs() for |Omega| and |Lambda|
    # min() for the minimum of J and (|Lambda| + S)
    rhs = min(J, abs(Lambda) + S)
    
    # Check the condition
    return abs(Omega) <= rhs
end
#
function generate_allowed_AM_values(state::Int64, J::Float64)
    S = (Potential[state].mult - 1)/2 # This correctly derives S from multiplicity
    #
    Sigma = (S == 0.0) ? [0.0] : [s for s=-S:S] # This correctly generates all valid Sigma projections for a given S
    #
    abs_Lambda_el = Potential[state].lambda 
    Lambda_vals_for_state = (abs_Lambda_el == 0.0) ? [0.0] : [-abs_Lambda_el, abs_Lambda_el]
    #
    Sig = [] # Initialize empty lists to store *selected* values
    #
    Lam = [] # Initialize empty lists to store *selected* values
    #
    Ome = [] # Initialize empty lists to store *selected* values
    #
    @views for L in Lambda_vals_for_state # Iterate through the *signed* Lambda values
        for ∑ in Sigma # Iterate through the Sigma values
            O = L+∑    # Calculate Omega
            # The check_omega_condition function should use the S_el (S for the current electronic state)
            if check_omega_condition(J, O, L, S) # Corrected arguments for check_omega_condition
                push!(Lam,L)
                push!(Sig,∑)
                push!(Ome,O)
            end
        end
    end
    #
    return Lam, Sig, Ome, S
end
#
function build_rotational_subBlocks(J)
    #
    ## initialise the number of electronic states
    Nstates = length( Calculation["method"].states )
    #
    ## initialise the list to hold rotational matrices
    Hrot_list = Vector{Matrix{Float64}}()
    for i=1:Nstates
        #
        Lambda, Sigma, Omega, S = generate_allowed_AM_values(i, J)
        #
        ## initialise the rotational kinetic energy matrix for the given electronic state 
        Nbasis = length(Omega)  # rotational subspace dimension for this state
        Hrot = zeros(Float64, Nbasis, Nbasis)  # empty matrix to fill in
        #
        ## compute J^2
        Jsqrd = J * (J + 1)
        #
        ## compute S^2
        Ssqrd = S * (S + 1)
        #
        for idx in eachindex(Lambda)
            for jdx=idx:lastindex(Lambda)
                Lambda_i, Sigma_i, Omega_i = Lambda[idx], Sigma[idx], Omega[idx]
                #
                Lambda_j, Sigma_j, Omega_j = Lambda[jdx], Sigma[jdx], Omega[jdx]
                #
                ## compute matrix elements in |J, Omega> basis
                #
                ## <J,Omega|Jz|J,Omega> = Omega
                Jz = (Omega_i == Omega_j) ? Omega_i : 0.0
                #
                ## < J, Omega -+ 1 | J± | J, Omega > = sqrt(J(J+1) - Omega(Omega -+ 1))
                if Omega_i == Omega_j - 1
                    Jmp = sqrt(Jsqrd - Omega_j*(Omega_j - 1))
                elseif Omega_i == Omega_j + 1
                    Jmp = sqrt(Jsqrd - Omega_j*(Omega_j + 1))
                else
                    Jmp = 0.0
                end
                #
                ## compute matrix elements in |Lambda, S, Sigma>  spin basis
                #
                ## <Lambda,S,∑|Sz|Lambda,S,∑> = ∑
                Sz = (Sigma_i == Sigma_j) ? Sigma_i : 0.0
                #
                ## < Lambda, S, ∑ ± 1 | S± | Lambda, S, ∑ > = sqrt(S(S+1) - Sigma(Sigma ± 1))
                if Sigma_i == Sigma_j + 1
                    Spm = sqrt(Ssqrd - Sigma_j*(Sigma_j + 1))
                elseif Sigma_i == Sigma_j - 1
                    Spm = sqrt(Ssqrd - Sigma_j*(Sigma_j - 1))
                else
                    Spm = 0.0
                end
                #
                ## populate Hrot matrix for the current electronic state
                diagonal_matel = (Jsqrd - Jz^2) + (Ssqrd - Sz^2)
                off_diagonal_matel = Jmp * Spm
                #
                matel = ((idx == jdx) ? diagonal_matel : 0.0) + off_diagonal_matel 
                #
                Hrot[idx,jdx] = matel
                Hrot[jdx,idx] = matel
            end
        end
        #
        push!(Hrot_list, Hrot)
    end
    #
    return Hrot_list
end
#
function generate_allowed_AM_values_2(state::Int64, J::Float64)
    S = (Potential[state].mult - 1)/2 # This correctly derives S from multiplicity
    #
    Sigma = (S == 0.0) ? [0.0] : [s for s=-S:S] # This correctly generates all valid Sigma projections for a given S
    #
    L = abs(Potential[state].lambda)
    #
    ## compute inversion QN
    if (L == 0)&(Potential[state].symmetry == "-")
        inv=1
    else
        inv=0
    end
    #
    Sig = [] # Initialize empty lists to store Sigma values
    #
    Lam = [] # Initialize empty lists to store Lambda values
    #
    Ome = [] # Initialize empty lists to store Omega values
    #
    phase = [] # Initialize empty lists to store phase values
    #
    @views for ∑ in Sigma # Iterate through the Sigma values
        O = L+∑           # Calculate Omega
        # The check_omega_condition function should use the S_el (S for the current electronic state)
        if check_omega_condition(J, O, L, S) # Corrected arguments for check_omega_condition
            #
            ## Compute ε phase
            exponent = Int(round(inv - L + S - ∑ + J - O))
            ε = (-1)^exponent
            #
            push!(Lam,L)
            push!(Sig,∑)
            push!(Ome,O)
            push!(phase,ε)
        end
    end
    #
    return Lam, Sig, Ome, S, phase
end
#
function build_rotational_parity_subBlocks(J, Tau)
    #
    ## initialise the number of electronic states
    Nstates = length( Calculation["method"].states )
    #
    ## initialise the list to hold rotational matrices
    Hrot_list = Vector{Matrix{Float64}}()
    for i=1:Nstates
        #
        Lambda, Sigma, Omega, S, phases = generate_allowed_AM_values_2(i, J)
        #
        ## initialise the rotational kinetic energy matrix for the given electronic state 
        Nbasis = length(Omega)  # rotational subspace dimension for this state
        Hrot = zeros(Float64, Nbasis, Nbasis)  # empty matrix to fill in
        #
        ## compute J^2
        Jsqrd = J * (J + 1)
        #
        ## compute S^2
        Ssqrd = S * (S + 1)
        #
        ## check if it is a singlet sigma + state
        # ...
        #
        ## compute rot matrix elements
        for idx in eachindex(Lambda)
            for jdx=idx:lastindex(Lambda) 
                Lambda_i, Sigma_i, Omega_i, phase_i = Lambda[idx], Sigma[idx], Omega[idx], phases[idx]
                Lambda_j, Sigma_j, Omega_j, phase_j = Lambda[jdx], Sigma[jdx], Omega[jdx], phases[jdx]
                #
                ## compute matrix elements in |Lambda, S, Sigma>  spin basis
                #
                ## <Lambda,S,∑|Sz^2|Lambda,S,∑> = ∑^2
                Sz2 = 0.0 #(Sigma_i == Sigma_j) ? Sigma_i^2 : 0.0
                #
                if Sigma_i == Sigma_j
                    Sz2 += Sigma_i^2
                elseif -Sigma_i == Sigma_j
                    Sz2 += Tau * Sigma_i^2
                end
                #
                ## compute matrix elements in |J, Omega> basis
                #
                ## <J,Omega|Jz|J,Omega> = Omega = 0 due to symmeterisation
                Jz2 = (Omega_i == Omega_j) ? Omega_i^2 : 0.0
                #
                ## compute spin-uncoupling ladder operators
                #
                ## < Lambda, S, ∑ ± 1 | S± | Lambda, S, ∑ > = sqrt(S(S+1) - Sigma(Sigma ± 1))
                Sm = 0.0
                Sp = 0.0
                #
                ## < J, Omega -+ 1 | J± | J, Omega > = sqrt(J(J+1) - Omega(Omega -+ 1))
                Jm = 0.0
                Jp = 0.0
                #
                JmSp = 0.0
                JpSm = 0.0
                #
                ## compute J+S-
                println(Omega_i," ",Omega_j)
                if (Omega_i == Omega_j - 1)&(Sigma_i == Sigma_j - 1)
                    Jp += 0.5 * sqrt(Jsqrd - Omega_j*(Omega_j - 1))
                    Sm += 0.5 * sqrt(Ssqrd - Sigma_j*(Sigma_j - 1))
                    #
                    JpSm += Jp * Sm
                end
                #
                if (Omega_i == -Omega_j - 1)&(Sigma_i == -Sigma_j - 1)
                    Jp += 0.5 * Tau * sqrt(Jsqrd - (-Omega_j)*((-Omega_j) - 1))
                    Sm += 0.5 * Tau * sqrt(Ssqrd - (-Sigma_j)*((-Sigma_j) - 1))
                    #
                    JpSm += Jp * Sm
                end
                #
                if (-Omega_i == Omega_j - 1)&(-Sigma_i == Sigma_j - 1)
                    Jp += 0.5 * Tau * sqrt(Jsqrd - (Omega_j)*((Omega_j) - 1))
                    Sm += 0.5 * Tau * sqrt(Ssqrd - (Sigma_j)*((Sigma_j) - 1))
                    #
                    JpSm += Jp * Sm
                end
                #
                if (-Omega_i == -Omega_j - 1)&(-Sigma_i == -Sigma_j - 1)
                    Jp += 0.5 * sqrt(Jsqrd - (-Omega_j)*((-Omega_j) - 1))
                    Sm += 0.5 * sqrt(Ssqrd - (-Sigma_j)*((-Sigma_j) - 1))
                    #
                    JpSm += Jp * Sm
                end
                #
                ## compute J-S+
                if (Omega_i == Omega_j + 1)&(Sigma_i == Sigma_j + 1)
                    Jm += 0.5 * sqrt(Jsqrd - Omega_j*(Omega_j + 1))
                    Sp += 0.5 * sqrt(Ssqrd - Sigma_j*(Sigma_j + 1))
                    #
                    JmSp += Jm * Sp
                end
                #
                if (Omega_i == -Omega_j + 1)&(Sigma_i == -Sigma_j + 1)
                    Jm += 0.5 * Tau * sqrt(Jsqrd - (-Omega_j)*((-Omega_j) + 1))
                    Sp += 0.5 * Tau * sqrt(Ssqrd - (-Sigma_j)*((-Sigma_j) + 1))
                    #
                    JmSp += Jm * Sp
                end
                #
                if (-Omega_i == Omega_j + 1)&(-Sigma_i == Sigma_j + 1)
                    Jm += 0.5 * Tau * sqrt(Jsqrd - (Omega_j)*((Omega_j) + 1))
                    Sp += 0.5 * Tau * sqrt(Ssqrd - (Sigma_j)*((Sigma_j) + 1))
                    #
                    JmSp += Jm * Sp
                end
                #
                if (-Omega_i == -Omega_j + 1)&(-Sigma_i == -Sigma_j + 1)
                    Jm += 0.5 * sqrt(Jsqrd - (-Omega_j)*((-Omega_j) + 1))
                    Sp += 0.5 * sqrt(Ssqrd - (-Sigma_j)*((-Sigma_j) + 1))
                    #
                    JmSp += Jm * Sp
                end                
                # #
                # ## < J, Omega -+ 1 | J± | J, Omega > = sqrt(J(J+1) - Omega(Omega -+ 1))
                # Jm = 0.0
                # Jp = 0.0
                # #
                # ## hardcode elements : J+
                # println(Omega_i," ",Omega_j)
                # if Omega_i == Omega_j - 1
                #     Jp += 0.5 * sqrt(Jsqrd - Omega_j*(Omega_j - 1))
                # end
                # # println(Omega_i == Omega_j - 1, " ",Jp," ",0.5 * Tau * sqrt(Jsqrd - (Omega_j)*((Omega_j) - 1)))
                # println(Jp)
                # #
                # if -Omega_i == Omega_j - 1
                #     Jp += 0.5 * Tau * sqrt(Jsqrd - Omega_j*(Omega_j - 1))
                # end
                # println(Jp)
                # #
                # if Omega_i == -Omega_j - 1
                #     Jp += 0.5 * Tau * sqrt(Jsqrd - (-Omega_j)*((-Omega_j) - 1))
                # end
                # println(Jp)
                # #
                # if -Omega_i == -Omega_j - 1
                #     Jp += 0.5  * sqrt(Jsqrd - (-Omega_j)*((-Omega_j) - 1))
                # end
                # println(Jp)
                # #
                # ## hardcode elements : J-
                # if Omega_i == Omega_j + 1
                #     Jm += 0.5 * sqrt(Jsqrd - Omega_j*(Omega_j + 1))
                # end
                # #
                # if -Omega_i == Omega_j + 1
                #     Jm += 0.5 * Tau * sqrt(Jsqrd - Omega_j*(Omega_j + 1))
                # end
                # #
                # if Omega_i == -Omega_j + 1
                #     Jm += 0.5 * Tau * sqrt(Jsqrd - (-Omega_j)*((-Omega_j) + 1))
                # end
                # #
                # if -Omega_i == -Omega_j + 1
                #     Jm += 0.5 * sqrt(Jsqrd - (-Omega_j)*((-Omega_j) + 1))
                # end
                # #
                # println("Jp = ",Jp, " Jm = ",Jm)
                # #
                # ## < Lambda, S, ∑ ± 1 | S± | Lambda, S, ∑ > = sqrt(S(S+1) - Sigma(Sigma ± 1))
                # Sm = 0.0
                # Sp = 0.0
                # #
                # ## hardcode elements : S+
                # if Sigma_i == Sigma_j + 1
                #     Sp += 0.5 * sqrt(Ssqrd - Sigma_j*(Sigma_j + 1))
                # end
                # #
                # if -Sigma_i == Sigma_j + 1
                #     Sp += 0.5 * Tau * sqrt(Ssqrd - Sigma_j*(Sigma_j + 1))
                # end
                # #
                # if Sigma_i == -Sigma_j + 1
                #     Sp += 0.5 * Tau * sqrt(Ssqrd - (-Sigma_j)*((-Sigma_j) + 1))
                # end
                # #
                # if -Sigma_i == -Sigma_j + 1
                #     Sp += 0.5 * sqrt(Ssqrd - (-Sigma_j)*((-Sigma_j) + 1))
                # end
                # #
                # ## hardcode elements : S-
                # if Sigma_i == Sigma_j - 1
                #     Sm += 0.5 * sqrt(Ssqrd - Sigma_j*(Sigma_j - 1))
                # end
                # #
                # if -Sigma_i == Sigma_j - 1
                #     Sm += 0.5 * Tau * sqrt(Ssqrd - Sigma_j*(Sigma_j - 1))
                # end
                # #
                # if Sigma_i == -Sigma_j - 1
                #     Sm += 0.5 * Tau * sqrt(Ssqrd - (-Sigma_j)*((-Sigma_j) - 1))
                # end
                # #
                # if -Sigma_i == -Sigma_j - 1
                #     Sm += 0.5 * sqrt(Ssqrd - (-Sigma_j)*((-Sigma_j) - 1))
                # end
                #
                ## Hrot = (J^2 - Jz) + (S^2 -Sz) - (J+S- +J-S+) + O(L)
                #
                ## populate Hrot matrix for the current electronic state
                diagonal_matel = (Jsqrd - Jz2) + (Ssqrd - Sz2)
                off_diagonal_matel = - (JpSm + JmSp)
                println(off_diagonal_matel)
                println()
                #
                matel = ((idx == jdx) ? diagonal_matel : 0.0) + off_diagonal_matel 
                #
                Hrot[idx,jdx] = matel
                Hrot[jdx,idx] = matel
            end
        end
        #
        push!(Hrot_list, Hrot)
    end
    #
    return Hrot_list
end
#
function build_Hrot(B::Matrix{Float64}, Hrot_list::Vector{Matrix{Float64}})
    #
    ## initialise contracted vibronic basis size
    vmax = Calculation["vibronic_solver"].contraction
    #
    ## initialise rotational block sizes
    rot_dims = [size(Hrot_list[i], 1) for i in 1:length(Hrot_list)]
    #
    ## initialise total rovibronic Hamiltonian dimension
    tot_rovibronic_dim = sum([vmax[state]*rot_dims[state] for state=1:lastindex(vmax)])
    #
    ## initialise rovibronic Hamiltonian
    Hrot = zeros(Float64, tot_rovibronic_dim, tot_rovibronic_dim)
    #
    ##
    Nstates = length( Calculation["method"].states )
    for state=1:Nstates
        #
        ## prepare vibrational block of <1/r^2>
        i = serialise_electronic_vib_indices(state, 1, vmax)
        j = serialise_electronic_vib_indices(state, vmax[state], vmax)
        #
        ## extract <state, v| B | state'', v''> block
        B_state = B[i:j, i:j]
        #
        ## tensor product the B_block with the rotational hamiltonian in the rot basis
        Hrot_state_block = kron(B_state, Hrot_list[state])
        #
        ## compute the final rovibronic index
        start_idx = serialise_vibronic_rotational_indices(state, 1, 1, vmax, rot_dims)
        end_idx   = serialise_vibronic_rotational_indices(state, vmax[state], rot_dims[state], vmax, rot_dims)
        #
        ## populate the rovibronic hamiltonian in the full basis
        Hrot[start_idx:end_idx, start_idx:end_idx] = Hrot_state_block
    end
    #
    ##
    kefac = KE_factor(Calculation["method"].atoms...)
    return  Hrot * kefac, rot_dims, tot_rovibronic_dim
end
#
function build_Hvibronic_embedding(T::Matrix{Float64}, V::Matrix{Float64}, spin_rot_dims::Vector{Int64}, tot_rovibronic_dim::Int64)::Matrix{Float64}
    #
    ## extract block from T + V
    function extract_T_V_block!(Hvib, bra_state, ket_state, vmax, spin_rot_dims, Irot, Hrot)
            #
            ## prepare state vibrational block
            bra_i = serialise_electronic_vib_indices(bra_state, 1, vmax)
            bra_j = serialise_electronic_vib_indices(bra_state, vmax[bra_state], vmax)
            #
            ket_i = serialise_electronic_vib_indices(ket_state, 1, vmax)
            ket_j = serialise_electronic_vib_indices(ket_state, vmax[ket_state], vmax)
            #
            ## extract <state, v| T + V | state'', v''> block
            T_V_block = Hvib[bra_i:bra_j, ket_i:ket_j]
            #
            ## tensor product the T_V_block with the rotational hamiltonian in the rot basis
            Hrovib_state_block = kron(T_V_block, Irot)
            #
            ## compute the final rovibronic index
            bra_start_idx = serialise_vibronic_rotational_indices(bra_state, 1, 1, vmax, spin_rot_dims)
            bra_end_idx   = serialise_vibronic_rotational_indices(bra_state, vmax[bra_state], spin_rot_dims[bra_state], vmax, spin_rot_dims)
            #
            ket_start_idx = serialise_vibronic_rotational_indices(ket_state, 1, 1, vmax, spin_rot_dims)
            ket_end_idx   = serialise_vibronic_rotational_indices(ket_state, vmax[ket_state], spin_rot_dims[ket_state], vmax, spin_rot_dims)
            #
            ## populate the rovibronic hamiltonian in the full basis
            Hrot[bra_start_idx:bra_end_idx, ket_start_idx:ket_end_idx] = Hrovib_state_block
    end
    #
    ## initialise the vmax embedding
    vmax = Calculation["vibronic_solver"].contraction
    #
    ## initialise total rovibronic Hamiltonian
    Hvibronic_embedding = zeros(Float64, tot_rovibronic_dim, tot_rovibronic_dim)
    #
    ## initialise the total vibronic Hamiltonian
    Hvib = T + V
    #
    ## prepare state vibrational block and map it to the full rovibronic Hamiltonian
    Nstates = length( Calculation["method"].states )
    for bra_state=1:Nstates
        for ket_state=bra_state:Nstates
            #
            ## Only process if symmetry and calculation conditions are satisfied
            if !check_symmetry_and_states(bra_state, ket_state)
                continue
            end
            #
            ## compute the spin-rotational identity
            Irot = Matrix{Float64}(I, spin_rot_dims[bra_state], spin_rot_dims[ket_state])
            #
            ## upper triangle
            extract_T_V_block!(Hvib, bra_state, ket_state, vmax, spin_rot_dims, Irot, Hvibronic_embedding) 
            #
            ## lower triangle
            if bra_state ≠ ket_state
                extract_T_V_block!(Hvib, ket_state, bra_state, vmax, spin_rot_dims, Irot, Hvibronic_embedding)
            end
        end
    end
    #
    return Hvibronic_embedding
end
#
function quantum_number_bookkeeping(J::Float64, Tau::Int)::Tuple{Dict{Int64, Vector{Number}},Dict{String , Int}}
    #
    ## electronic states
    el_states = Calculation["method"].states
    #
    ## vibrational index's
    vmax = Calculation["vibronic_solver"].contraction
    #
    ## spin-rotational
    spin_rot_quanta = [generate_allowed_AM_values_2(state, J) for state in el_states]
    spin_rot_dims = [length(s[1]) for s in spin_rot_quanta]
    #
    ## initialise bookkeeping dictionary
    QN_book = Dict()
    #
    ##
    for state in el_states
        #
        Lambda, Sigma, Omega, S, phase = spin_rot_quanta[state]
        #
        for v_idx=1:vmax[state]
            #
            for rot_idx=1:spin_rot_dims[state]
                global_index = serialise_vibronic_rotational_indices(state, v_idx, rot_idx, vmax, spin_rot_dims)
                #
                QN_book[global_index] = [state,
                                         v_idx-1,
                                         Lambda[rot_idx],
                                         Sigma[rot_idx],
                                         J,
                                         Omega[rot_idx],
                                         Tau]
            end
        end
    end
    #
    ## create dictionary for QN positions in array
    quantum_number_map = Dict("state"=>1,
                  "v"=>2,
                  "Lambda"=>3,
                  "Sigma"=>4,
                  "J"=>5,
                  "Omega"=>6,
                  "Parity"=>7)
    #
    return QN_book, quantum_number_map
end

#
function quantum_number_bookkeeping_old(J::Float64)::Tuple{Dict{Int64, Vector{Number}},Dict{String , Int64}}
    #
    ## electronic states
    el_states = Calculation["method"].states
    #
    ## vibrational index's
    vmax = Calculation["vibronic_solver"].contraction
    #
    ## spin-rotational
    spin_rot_quanta = [generate_allowed_AM_values(state, J) for state in el_states]
    spin_rot_dims = [length(s[1]) for s in spin_rot_quanta]
    #
    ## initialise bookkeeping dictionary
    QN_book = Dict()
    #
    ##
    for state in el_states
        #
        Lambda, Sigma, Omega, S = spin_rot_quanta[state]
        #
        S = (Potential[state].mult - 1)/2
        #
        for v_idx=1:vmax[state]
            #
            for rot_idx=1:spin_rot_dims[state]
                global_index = serialise_vibronic_rotational_indices(state, v_idx, rot_idx, vmax, spin_rot_dims)
                #
                QN_book[global_index] = [state,
                                         v_idx-1,
                                         Lambda[rot_idx],
                                         Sigma[rot_idx],
                                         J,
                                         Omega[rot_idx]]
            end
        end
    end
    #
    ## create dictionary for QN positions in array
    quantum_number_map = Dict("state"=>1,
                  "v"=>2,
                  "Lambda"=>3,
                  "Sigma"=>4,
                  "J"=>5,
                  "Omega"=>6)
    #
    return QN_book, quantum_number_map
end


# function build_rotational_subBlocks(J)
#     #
#     ## initialise the number of electronic states
#     Nstates = length( Calculation["method"].states )
#     #
#     ## initialise the list to hold rotational matrices
#     Hrot_list = Vector{Matrix{Float64}}()
#     for i=1:Nstates
#         for j=i:Nstates
#             #
#             iLambda, iSigma, iOmega, iS = generate_allowed_AM_values(i, J)
#             jLambda, jSigma, jOmega, jS = generate_allowed_AM_values(i, J)
#             #
#             ## initialise the rotational kinetic energy matrix for the given electronic state 
#             Nbasis_i = length(iOmega)  # rotational subspace dimension for this state
#             Nbasis_j = length(jOmega)  # rotational subspace dimension for this state
#             #
#             Hrot = zeros(Float64, Nbasis_i, Nbasis_j)  # empty matrix to fill in
#             #
#             ## if spin-rotationally diagonal
#             if (iLambda, iSigma)
#             #
#             ## compute J^2
#             Jsqrd = J * (J + 1)
#             #
#             ## compute S^2
#             Ssqrd = S * (S + 1)
#             #
#             for idx in eachindex(Lambda)
#                 for jdx=idx:lastindex(Lambda)
#                     Lambda_i, Sigma_i, Omega_i = Lambda[idx], Sigma[idx], Omega[idx]
#                     #
#                     Lambda_j, Sigma_j, Omega_j = Lambda[jdx], Sigma[jdx], Omega[jdx]
#                     #
#                     ## compute matrix elements in |J, Omega> basis
#                     #
#                     ## <J,Omega|Jz|J,Omega> = Omega
#                     Jz = (Omega_i == Omega_j) ? Omega_i : 0.0
#                     #
#                     ## < J, Omega -+ 1 | J± | J, Omega > = sqrt(J(J+1) - Omega(Omega -+ 1))
#                     if Omega_i == Omega_j - 1
#                         Jmp = sqrt(Jsqrd - Omega_j*(Omega_j - 1))
#                     elseif Omega_i == Omega_j + 1
#                         Jmp = sqrt(Jsqrd - Omega_j*(Omega_j + 1))
#                     else
#                         Jmp = 0.0
#                     end
#                     #
#                     ## compute matrix elements in |Lambda, S, Sigma>  spin basis
#                     #
#                     ## <Lambda,S,∑|Sz|Lambda,S,∑> = ∑
#                     Sz = (Sigma_i == Sigma_j) ? Sigma_i : 0.0
#                     #
#                     ## < Lambda, S, ∑ ± 1 | S± | Lambda, S, ∑ > = sqrt(S(S+1) - Sigma(Sigma ± 1))
#                     if Sigma_i == Sigma_j + 1
#                         Spm = sqrt(Ssqrd - Sigma_j*(Sigma_j + 1))
#                     elseif Sigma_i == Sigma_j - 1
#                         Spm = sqrt(Ssqrd - Sigma_j*(Sigma_j - 1))
#                     else
#                         Spm = 0.0
#                     end
#                     #
#                     ## populate Hrot matrix for the current electronic state
#                     diagonal_matel = (Jsqrd - Jz^2) + (Ssqrd - Sz^2)
#                     off_diagonal_matel = Jmp * Spm
#                     #
#                     matel = ((idx == jdx) ? diagonal_matel : 0.0) + off_diagonal_matel 
#                     #
#                     Hrot[idx,jdx] = matel
#                     Hrot[jdx,idx] = matel
#                 end
#             end
#             #
#             push!(Hrot_list, Hrot)
#         end
#     end
#     #
#     return Hrot_list
# end