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
function Built_nonadiabatic_T_V() #:: TYPE DECLARE NEEDED
    #
    ## flatten multidimensional electronic+vibrational indexing
    function serialise_electronic_vib_indices(state, v_idx, contraction_array)
        return ( state - 1 ) * contraction_array[ state ] + v_idx
    end
    #
    ## initialise contracted vibronic basis size
    vmax = Calculation["vibronic_solver"].contraction
    contracted_vibronic_dim = sum(vmax)
    #
    ## initialise the kinetic energy and electronic Hamiltonians
    T = zeros(Float64, contracted_vibronic_dim, contracted_vibronic_dim)
    V = zeros(Float64, contracted_vibronic_dim, contracted_vibronic_dim)
    #
    ## initialise the number of electronic states
    Nstates = length( Calculation["method"].states )
    #
    ## initialise the NAC squared matrix: second DDR
    W2 = zeros(Float64, length(r), Nstates, Nstates)
    K = map(idx -> NACMat[idx,:,:] * NACMat[idx,:,:], collect(1:lastindex(r)))
    for i=1:Nstates
        for j=i:Nstates
            W2[:,i,j] .= [K[idx][i,j] for idx=1:lastindex(r)]
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
        for v_idx=1:vmax[i]
            for v_jdx=v_idx:vmax[i]
                #
                row    = serialise_electronic_vib_indices(i, v_idx, vmax)
                column = serialise_electronic_vib_indices(i, v_jdx, vmax)
                #
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
    return T, V
end
#
