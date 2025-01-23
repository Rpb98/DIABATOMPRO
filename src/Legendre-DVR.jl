using LegendrePolynomials
using PyPlot
using FastGaussQuadrature, LinearAlgebra
using LinearAlgebra

# Example: Compute roots and weights
r = LinRange(-1,1,1000)
N = 5
Pn = Pl.(r,100)
x, w = gausslegendre(100)

# f(x) = x^4
# I = dot(w, f.(x))
# I ≈ 2/5
# roots = Polynomials.roots(P)
# weights = legendre_weights(roots, N)
# println("Roots: ", roots)
# println("Weights: ", weights)

# plt.figure()
# plt.plot(r,Pn)
#
function HarmonicOscillator(r; re = 10.5, k=1000)
    return 0.5*k*(r-re)^2
end

R = LinRange(0.5,20.5,1000)

plt.figure()
plt.plot(R,HarmonicOscillator.(R,re=10.5,k=1000))
#
function quad2bond(xi, R_min, R_max)
    # Scale the Legendre roots from [-1, 1] to [R_min, R_max]
    return 0.5 * ((R_max - R_min) * xi .+ (R_max + R_min))
end
#
function DVR_KinMatel(v, root_i, root_j)
    #
    ## compute roots and quadrature weights for legendre polynomial v 
    xv, wv = gausslegendre(v)       
    #
    ## check diagonal or off-diagonal DVR matrix element 
    if root_i != root_j                                                         # off-diagonal elements
        #
        xi, xj = xv[root_i], xv[root_j] # ith root and weight
        wi, wj = wv[root_i], wv[root_j] # jth root and weight
        #
        return sqrt(wi*wj)/(xi-xj)^2
    elseif root_i == root_j                                                     # diagonal elements
        #
        Nroots = length(xv)
        #
        Tii = 0
        #
        for i=1:Nroots
            for j=1:Nroots
                if i != j
                    #
                    xi, xj = xv[i], xv[j] # ith root and weight
                    wi, wj = wv[i], wv[j] # jth root and weight
                    #
                    Tii += sqrt(wi*wj)/(xi-xj)^2
                end
            end
        end
        #
        return Tii
    end 
end
#
function initialise_DVR_matrix(vmax)
    #
    ## compute number of roots of the vmax vibrational basis functions
    #
    N = 1
    #
    if vmax == 1
        N = 1
    else 
        N = Int(vmax*(vmax+1)/2)
    end
    #
    return zeros(Float64, N,N)
end
#
function build_DVR_KinMat(vmax)
    #
    ## initialise the matrix
    T = initialise_DVR_matrix(vmax)
    #
    ## populate the matrix for vibrational states up to vmax
    idx_shift = 0
    for v=1:vmax
        #
        ## number of DVR roots equals the vibrational QN
        Nv = v 
        #
        idx_shift += (v-1)
        #
        for i=1:Nv
            for j=i:Nv
                Tij = DVR_KinMatel(v, i, j)
                #
                ## Fill the matrix symmetrically
                T[i + idx_shift, j + idx_shift] = Tij
                if i != j
                    T[j + idx_shift, i + idx_shift] = Tij  # Symmetric entry for Tij
                end
            end
        end
    end
    #
    return T
end
#
function build_DVR_PotMat(vmax, fPot, r_min, r_max)
    #
    ## initialise the matrix
    V = initialise_DVR_matrix(vmax)
    #
    ## populate the matrix for vibrational states up to vmax
    idx_shift = 0
    for v=1:vmax
        #
        ## number of DVR roots equals the vibrational QN
        Nv = v 
        #
        idx_shift += (v-1)
        #
        ## compute quadrature roots and weights
        xv, wv = gausslegendre(v)
        #
        for i=1:Nv
            #
            ## scale the quadrature root to the bond length region
            ri = quad2bond(xv[i], r_min, r_max)
            #
            ## evaluate the potential at the weight
            Vii = fPot(ri)
            #
            ## Fill the matrix symmetrically
            V[i + idx_shift, i + idx_shift] = Vii
        end
    end
    #
    return V
end

T = build_DVR_KinMat(10)
V = build_DVR_PotMat(10,HarmonicOscillator,0.5,20.5)


function solve_vibrational_schrodinger(T, V, r)
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
    ## compute the vibronic wavefunctions from the DVR coefficients in eigenvectors
    # wavefunctions = []
    # for v=1:Nv
    #     coeff = eigenvectors[:,v]
    #     #
    #     for c in coeff

    
    return eigenvalues, eigenvectors
end

function DVRtoBasis(r,v,root, r_min, r_max)
    xv, wv = gausslegendre(v)
    #
    ri = quad2bond(xv[root],r_min, r_max)
    #
    x = (2 * r - (r_min + r_max)) / (r_max - r_min)
    #
    if r == ri
        return 1
    else
        return Pl(x,v)/((r-ri)*dnPl(x,v,1))
    end
end

# Example usage:
# Assuming `T` and `V` are already computed:
E, Ψ = solve_vibrational_schrodinger(T, V, r)

# Print the vibrational energy levels
println("Vibrational energy levels: ", E)

plt.figure()
plt.plot(R, DVRtoBasis.(R,3,1,0.5,20.5))

# Ψ contains the wavefunctions in the DVR basis
