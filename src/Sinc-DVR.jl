using LegendrePolynomials
using PyPlot
using FastGaussQuadrature
using LinearAlgebra
using Dierckx

function KE_factor(m1,m2)
    #
    ## constants
    h = 6.626069570000*10^(-34)   # J s
    c = 299792458.0000       # m/s
    amu = 1.660538921000*10^(-27) # kg
    #
    ## reduced mass
    mu = amu*m1*m2/(m1+m2)
    #
    ## KE factor
    factor = (h/(8*pi^2*mu*c))*10^(18)
    #
    return factor
end

global state_v_wav

Ne = 2
vmax = 500
Ngrid = 10001
state_v_wav = Array{Float64}(undef, Ne, vmax, Ngrid)


# Example: Compute roots and weights
rmin = 0.5
rmax = 6
r = LinRange(rmin,rmax,Ngrid)
h = r[2]-r[1]
#
function sinc_function(x, xj, h)
    # Avoid division by zero: handle x == x_j separately
    if abs(x - xj) < 1e-10
        return 1.0
    else
        return sin(pi * (x - xj) / h) / (pi * (x - xj) / h)
    end
end
#
function HarmonicOscillator(r; re = 3.5, nu = 1350, m1 = 12, m2 = 1.007825032, Ve = 0)
    h = 6.626069570000*10^(-34)   # J s
    c = 299792458.0000            # m/s
    amu = 1.660538921000*10^(-27) # kg
    #
    ## reduced mass
    mu = amu*m1*m2/(m1+m2)
    #
    ##
    factor = 4 * pi^2 * c * 10^(-18) / h
    return Ve + factor * 0.5 * mu * nu^2 * (r-re)^2
end
#
function SincDVR_KinMat(vmax,r_min, r_max, factor)
    #
    T = zeros(Float64, vmax, vmax)
    #
    ## compute quadrature points
    R = LinRange(r_min, r_max,vmax)
    #
    h = R[2] - R[1]
    #
    for i=1:vmax
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
function SincDVR_PotMat(vmax, fPot, r_min, r_max)
    #
    ## evaluate quadrature points
    R = LinRange(r_min,r_max,vmax)
    #
    ## compute potential on the qwuadrature points
    V = Diagonal(fPot.(R))
    #
    return V
end
#
function solve_vibrational_schrodinger(T, V, r, R, h, electronic_state)
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
    ## precompute Sinc-DVR basis functions
    φ = [sinc_function.(r, R[i], h) for i=1:lastindex(R)]
    #
    ## compute uncoupled vibronic wavefunctions as basis and populate an array
    for idx=1:Nv
        #
        v = idx - 1
        #
        state_v_wav[electronic_state,v+1,:] .= Vibronic_Wavefunction(v, eigenvectors, φ, R, h, r)
    end
    #
    return eigenvalues, eigenvectors
end
#
function Vibronic_Wavefunction(v, eigenvectors, φ, Rquad, h, r; thresh = 1e-100)
    #
    ## sum the vibronic basis with the weights coeff, only compute if coefficient is above threshold
    Ψv = zeros(Float64, length(r))
    #
    coeff = eigenvectors[:,v+1]
    #
    for (i, c) in enumerate(coeff)
        Ψv .+= c .* φ[i]
    end
    #
    return Ψv
end
#
function adiabatic_HO(r, V1d, V2d, gamma, re)
    #
    ##
    beta = pi/4 .+ 0.5 .* atan.((r .- re) ./ gamma)
    W12 = 0.5 .* gamma .* ( (r.-re).^2 .+ gamma^2).^(-1)
    #
    ∂Vd = V2d .- V1d
    #
    DC = 0.5 .* tan.(2 .* beta) .* ∂Vd
    #
    ##
    V1a = (V1d .+ V2d) .- sqrt.(∂Vd.^2 + 4 .* DC.^2) #V1d .* cos.(beta).^2 .+ V2d .* sin.(beta).^2
    V2a = (V1d .+ V2d) .+ sqrt.(∂Vd.^2 + 4 .* DC.^2) #V2d .* cos.(beta).^2 .+ V1d .* sin.(beta).^2
    #
    return [V1a, V2a], W12, beta
end
#
vmax = 500
#
V1_ = HarmonicOscillator.(r)
V2_ = HarmonicOscillator.(r,re=4.2,nu=700,Ve=1000)
V__ = [V1_,V2_]

spl = Spline1D(r, V1_ .- V2_)
x0 = roots(spl)
rcross = x0[findfirst(3.5 .< x0 .< 4.5)]

#
## evaluate quadrature points
R = LinRange(r[1],r[end], vmax)
h = R[2] - R[1]
#
V1 = Diagonal(HarmonicOscillator.(R))
V2 = Diagonal(HarmonicOscillator.(R,re=4.2,nu=700,Ve=1000))
V_ = [V1,V2]

V, _, _ = adiabatic_HO(R, HarmonicOscillator.(R), HarmonicOscillator.(R,re=4.2,nu=700,Ve=1000), 0.03, rcross)

V = [Diagonal(V[1]), Diagonal(V[2])]

Va, W12, b12 = adiabatic_HO(r, V1_, V2_, 0.03, rcross)

T = SincDVR_KinMat(vmax, r[1], r[end], KE_factor(12,1.007825032))

plt.figure()
for (idx,Vi) in enumerate(V)
    E, eigenvectors = solve_vibrational_schrodinger(T, Vi, r, R, h, idx)
    
    plt.plot(r,Va[idx])
    for (v,e) in enumerate(E[1:10])
        Psi_v = state_v_wav[idx,v,:] #Vibronic_Wavefunction(v-1, eigenvectors, r[1], r[end], vmax, r)
        # plt.axhline(e)
        #
        ## compute the scaling for the wavefunction to plot
        ∂E = E[v+1] - e
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
        plt.plot(r[idx_i:idx_f], scale .* prob_density_norm[idx_i:idx_f] .+ e)
    end
    plt.ylim(0,E[5]*1.1)
end

# # #
# # ## now compute the adiabatic Hamiltonian in vibronic basis
# # Hnac = Array{Float64}(undef, Ne, Ne, Ngrid)
# k1 = 50000
# k2 = 40000
# r01 = 2.5
# r02 = 4.5
# #
# # k2 = 0.5*k1*r01/r02
# # #
# # f = -(k2*r02)/(k1*r01 - k2*r02)
# #
# V1 = 0.5*k1 .* (r .- r01).^2
# V2 = 0.5*k2 .* (r .- r02).^2
# #
# A = 0.5*k1
# B = 0.5*k2



# var =(2*A .* (r .- r01).^2) ./ (B .* (r .- r02).^2) .- 1

# # alpha = log.(var) ./ (r .- 3.4)


# # V = f .* V2 .+ (1-f) .* V1
# alpha = -10
# rc = 3.5

# gamma = acosh.(sqrt.(cosh.(r .- r02) .* (1 .+ exp.(alpha .* (r .- rc)))))

# f = cosh.(gamma).^(-2)
# # k2 = LinRange(0.1,100000,100000)
# # f = -(k2 .* r02)./(k1*r01 .- k2 .*r02)
# plt.figure()
# plt.plot(r,f)
# # plt.yscale("symlog")
# # plt.plot(k2,f)
# # plt.ylim(0,50000)

# # EHO = [(v-1+0.5)*1350 for v=1:vmax]

# # vs = [v for v=1:vmax]
# # plt.figure()
# # plt.plot(vs,abs.(E .- EHO))
# # plt.yscale("log")
# # plt.plot(vs,EHO)




# # plt.figure()
# # plt.plot([i-1 for i=1:vmax], E)



# # # R = LinRange(0.5,20.5,1000)

# # # plt.figure()
# # # plt.plot(R,HarmonicOscillator.(R,re=10.5,k=1000))
# # #
# # function quad2bond(xi, R_min, R_max)
# #     # Scale the Legendre roots from [-1, 1] to [R_min, R_max]
# #     return 0.5 * ((R_max - R_min) * xi .+ (R_max + R_min))
# # end
# # #
# # function DVR_KinMatel(v, root_i, root_j)
# #     #
# #     ## compute roots and quadrature weights for legendre polynomial v 
# #     xv, wv = gausslegendre(v)       
# #     #
# #     ## check diagonal or off-diagonal DVR matrix element 
# #     if root_i != root_j                                                         # off-diagonal elements
# #         #
# #         xi, xj = xv[root_i], xv[root_j] # ith root and weight
# #         wi, wj = wv[root_i], wv[root_j] # jth root and weight
# #         #
# #         return sqrt(wi*wj)/(xi-xj)^2
# #     elseif root_i == root_j                                                     # diagonal elements
# #         #
# #         Nroots = length(xv)
# #         #
# #         Tii = 0
# #         #
# #         for i=1:Nroots
# #             for j=1:Nroots
# #                 if i != j
# #                     #
# #                     xi, xj = xv[i], xv[j] # ith root and weight
# #                     wi, wj = wv[i], wv[j] # jth root and weight
# #                     #
# #                     Tii += sqrt(wi*wj)/(xi-xj)^2
# #                 end
# #             end
# #         end
# #         #
# #         return Tii
# #     end 
# # end
# # #
# # function initialise_DVR_matrix(vmax)
# #     #
# #     ## compute number of roots of the vmax vibrational basis functions
# #     #
# #     N = 1
# #     #
# #     if vmax == 1
# #         N = 1
# #     else 
# #         N = Int(vmax*(vmax+1)/2)
# #     end
# #     #
# #     return zeros(Float64, N,N)
# # end
# # #
# # function build_DVR_KinMat(vmax)
# #     #
# #     ## initialise the matrix
# #     T = initialise_DVR_matrix(vmax)
# #     #
# #     ## populate the matrix for vibrational states up to vmax
# #     idx_shift = 0
# #     for v=1:vmax
# #         #
# #         ## number of DVR roots equals the vibrational QN
# #         Nv = v 
# #         #
# #         idx_shift += (v-1)
# #         #
# #         for i=1:Nv
# #             for j=i:Nv
# #                 Tij = DVR_KinMatel(v, i, j)
# #                 #
# #                 ## Fill the matrix symmetrically
# #                 T[i + idx_shift, j + idx_shift] = Tij
# #                 if i != j
# #                     T[j + idx_shift, i + idx_shift] = Tij  # Symmetric entry for Tij
# #                 end
# #             end
# #         end
# #     end
# #     #
# #     return T
# # end
# # #
# # function build_DVR_PotMat(vmax, fPot, r_min, r_max)
# #     #
# #     ## initialise the matrix
# #     V = initialise_DVR_matrix(vmax)
# #     #
# #     ## populate the matrix for vibrational states up to vmax
# #     idx_shift = 0
# #     for v=1:vmax
# #         #
# #         ## number of DVR roots equals the vibrational QN
# #         Nv = v 
# #         #
# #         idx_shift += (v-1)
# #         #
# #         ## compute quadrature roots and weights
# #         xv, wv = gausslegendre(v)
# #         #
# #         for i=1:Nv
# #             #
# #             ## scale the quadrature root to the bond length region
# #             ri = quad2bond(xv[i], r_min, r_max)
# #             #
# #             ## evaluate the potential at the weight
# #             Vii = fPot(ri)
# #             #
# #             ## Fill the matrix symmetrically
# #             V[i + idx_shift, i + idx_shift] = Vii
# #         end
# #     end
# #     #
# #     return V
# # end

# # T = build_DVR_KinMat(10)
# # V = build_DVR_PotMat(10,HarmonicOscillator,0.5,20.5)


# # function solve_vibrational_schrodinger(T, V, r)
# #     """
# #     Solve the vibrational Schrödinger equation Hψ = Eψ using DVR.
    
# #     Parameters:
# #         T: Kinetic energy matrix (square matrix)
# #         V: Potential energy matrix (diagonal matrix, same size as T)
    
# #     Returns:
# #         E: Eigenvalues (vibrational energy levels)
# #         Ψ: Eigenvectors (vibrational wavefunctions in DVR basis)
# #     """
# #     #
# #     ## Construct the total Hamiltonian
# #     H = T + V
# #     #
# #     ## Diagonalize the Hamiltonian
# #     eigenvalues, eigenvectors = eigen(H)
# #     #
# #     ## number of eigenvalues
# #     Nv = length(eigenvalues)
# #     #
# #     ## compute the vibronic wavefunctions from the DVR coefficients in eigenvectors
# #     # wavefunctions = []
# #     # for v=1:Nv
# #     #     coeff = eigenvectors[:,v]
# #     #     #
# #     #     for c in coeff

    
# #     return eigenvalues, eigenvectors
# # end

# # function DVRtoBasis(r,v,root, r_min, r_max)
# #     xv, wv = gausslegendre(v)
# #     #
# #     ri = quad2bond(xv[root],r_min, r_max)
# #     #
# #     x = (2 * r - (r_min + r_max)) / (r_max - r_min)
# #     #
# #     if r == ri
# #         return 1
# #     else
# #         return Pl(x,v)/((r-ri)*dnPl(x,v,1))
# #     end
# # end

# # # Example usage:
# # # Assuming `T` and `V` are already computed:
# # E, Ψ = solve_vibrational_schrodinger(T, V, r)

# # # Print the vibrational energy levels
# # println("Vibrational energy levels: ", E)

# # plt.figure()
# # plt.plot(R, DVRtoBasis.(R,3,1,0.5,20.5))

# # # Ψ contains the wavefunctions in the DVR basis
