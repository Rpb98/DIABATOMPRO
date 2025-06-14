using LegendrePolynomials
using PyPlot
using FastGaussQuadrature
using LinearAlgebra
using Dierckx
using Statistics

function KE_factor(m1::Float64,m2::Float64)::Float64
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
    return (h/(8*pi^2*mu*c))*10^(18)
end
#
function sinc_function(x::Float64, xj::Float64, h::Float64)::Float64
    # Avoid division by zero: handle x == x_j separately
    if abs(x - xj) < 1e-10
        return 1.0
    else
        return sin(pi * (x - xj) / h) / (pi * (x - xj) / h)
    end
end
#
function SincDVR_KinMat(vmax::Int,r_min::Float64, r_max::Float64, factor::Float64)::Matrix{Float64}
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
function solve_vibrational_schrodinger( T                :: Matrix{Float64}, 
                                        V                :: Diagonal{Float64, Vector{Float64}}, 
                                        r                :: Vector{Float64}, 
                                        h                :: Float64, 
                                        electronic_state :: Int64
                                        )                :: Tuple{Vector{Float64}, Matrix{Float64}}
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
    φ = [sinc_function.(r, x, h) for x in r]
    #
    ## compute uncoupled vibronic wavefunctions as basis and populate an array
    for idx=1:Nv
        #
        v = idx - 1
        #
        state_v_wav[electronic_state,v+1,:] .= Vibronic_Wavefunction(v, eigenvectors, φ, r)
    end
    #
    return eigenvalues, eigenvectors
end
#
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
    for (i, c) in enumerate(coeff)
        Ψv .+= c .* φ[i]
    end
    #
    return Ψv
end
#
@views function simps(f::Function, x::AbstractRange) #A range with odd number
    h = step(x)
    I= h/3*(f(x[1])+2*sum(f,x[3:2:end-2])+4*sum(f,x[2:2:end-1])+f(x[end]))
    return I
end

# from the function and range definition
function simps(f::Function, a::Real, b::Real, n::Integer) #n as an even number
    return simps(f, range(a, b, length=n+1))
end

@views function simps(fx::AbstractVector, h::Real)
    if length(fx)%2==1
        I= h/3*(fx[1]+2*sum(fx[3:2:end-2])+4*sum(fx[2:2:end-1])+fx[end])
    else
        I=h/3*(fx[1]+2*sum(fx[3:2:end-5])+4*sum(fx[2:2:end-4])+fx[end-3])+
        (3*h/8)*(fx[end-3]+3fx[end-2]+3fx[end-1]+fx[end])
    end
    return Float64(I)
end
# from the function values and range
function simps(fx::AbstractVector, a::Real, b::Real)
    return simps(fx, (b-a)/(length(fx)-1))
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

function vibronic_eigensolver(r, V, vmax)
    #
    ## convert potential object into vector of diagonal objects
    






global state_v_wav
Ne = 2
vmax = 1000
Ngrid = vmax
state_v_wav = Array{Float64}(undef, Ne, vmax, vmax)
# Example: Compute roots and weights
rmin = 0.5
rmax = 10
r = collect(LinRange(rmin,rmax,Ngrid))
h = r[2]-r[1]
#
V1 = Diagonal(HarmonicOscillator.(r))
V2 = Diagonal(HarmonicOscillator.(r,re=4.2,nu=700,Ve=1000))
V_ = [V1,V2]
#
T = SincDVR_KinMat(vmax, r[1], r[end], KE_factor(12.0,1.007825032))
#
V1_ = HarmonicOscillator.(r)
V2_ = HarmonicOscillator.(r,re=4.2,nu=700,Ve=1000)
#
V_ = [V1_,V2_]
#
V = [V1,V2]
plt.figure()

for (idx,Vi) in enumerate(V)
    E, eigenvectors = solve_vibrational_schrodinger(T, Vi, r, h, idx)
    
    plt.plot(R,V_[idx])
    for (v,e) in enumerate(E[1:5])
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
    plt.ylim(0,1e5)
end












# function SincDVR_KinMat(vmax::Int,r_min::Float64, r_max::Float64, factor::Float64)::Matrix{Float64}
#     #
#     T = zeros(Float64, vmax, vmax)
#     #
#     ## compute quadrature points
#     R = LinRange(r_min, r_max,vmax)
#     #
#     h = R[2] - R[1]
#     #
#     for i=1:vmax
#         for j=i:vmax
#             if i == j
#                 T[i,i] = pi^2 / (3 * h^2)
#             else
#                 T[i,j] = 2 * ( (-1)^(i-j) / (i-j)^2 ) / h^2
#                 T[j,i] = 2 * ( (-1)^(i-j) / (i-j)^2 ) / h^2
#             end
#         end
#     end
#     #
#     return factor * T
# end
# #
# function SincDVR_PotMat(vmax::Int, fPot, r_min::Float64, r_max::Float64)::Diagonal(Float64,Vector{Float64})
#     #
#     ## evaluate quadrature points
#     R = LinRange(r_min,r_max,vmax)
#     #
#     ## compute potential on the quadrature points
#     V = Diagonal(fPot.(R))
#     #
#     return V
# end
# #
# function solve_vibrational_schrodinger(T::Matrix{Float64}, V::Diagonal{Float64,Vector{Float64}}, r::Vector{Float64}, h::Float64, electronic_state::Int64)::Tuple{Vector{Float64}, Matrix{Float64}}

#     """
#     Solve the vibrational Schrödinger equation Hψ = Eψ using DVR.
    
#     Parameters:
#         T: Kinetic energy matrix (square matrix)
#         V: Potential energy matrix (diagonal matrix, same size as T)
    
#     Returns:
#         E: Eigenvalues (vibrational energy levels)
#         Ψ: Eigenvectors (vibrational wavefunctions in DVR basis)
#     """
#     #
#     ## Construct the total Hamiltonian
#     H = T + V
#     #
#     ## Diagonalize the Hamiltonian
#     eigenvalues, eigenvectors = eigen(H)
#     #
#     ## number of eigenvalues
#     Nv = length(eigenvalues)
#     #
#     ## precompute Sinc-DVR basis functions
#     φ = [sinc_function.(r, R[i], h) for i=1:lastindex(r)]
#     #
#     ## compute uncoupled vibronic wavefunctions as basis and populate an array
#     for idx=1:Nv
#         #
#         v = idx - 1
#         #
#         state_v_wav[electronic_state,v+1,:] .= Vibronic_Wavefunction(v, eigenvectors, φ, r)
#     end
#     #
#     return eigenvalues, eigenvectors
# end
# #
# function Vibronic_Wavefunction(v::Int64, eigenvectors::Matrix{Float64}, φ::Vector{Float64}, r::Vector{Float64}; thresh = 1e-100)::Vector{Float64}
#     #
#     ## sum the vibronic basis with the weights coeff, only compute if coefficient is above threshold
#     Ψv = zeros(Float64, length(r))
#     #
#     coeff = eigenvectors[:,v+1]
#     #
#     for (i, c) in enumerate(coeff)
#         Ψv .+= c .* φ[i]
#     end
#     #
#     return Ψv
# end
#
## now couple the vibronic states with NACs
#
# function FiniteDifference(x::Vector{Float64}, y::Vector{Float64}, d_order::Int64)::Vector{Float64}
#     local dy
#     #
#     ## derivative parameters for forward and backwards differences
#     assymetrical_first_derivative_parameters = [[-1,       1],
#                                                 [-1.5,     2, -0.5], 
#                                                 [-11/6,    3, -3/2,  1/3],
#                                                 [-25/12,   4, -3,	4/3,  -0.25],
#                                                 [-137/60,  5, -5,    10/3, -5/4,1 /5],
#                                                 [-49/20,   6, -15/2, 20/3, -15/4, 6/5,  -1/6],
#                                                 [-363/140, 7, -21/2, 35/3, -35/4, 21/5, -7/6,  1/7],
#                                                 [-761/280, 8, -14,   56/3, -35/2, 56/5, -14/3, 8/7, -1/8]]
#     assymetrical_second_derivative_parameters = [[1,          -2,      1],	 	 	 	 	 	 
#                                                  [2,          -5,      4,      -1],	 	 	 	 	 
#                                                  [35/12,      -26/3,   19/2,   -14/3,    11/12], 	 	 	 
#                                                  [15/4,       -77/6,   107/6,  -13,      61/12, -5/6],	 	 	 
#                                                  [203/45,     -87/5,   117/4,  -254/9,   33/2,  -27/5,   137/180],	 	 
#                                                  [469/90,     -223/10, 879/20, -949/18,  41,    -201/10, 1019/180, -7/10],
#                                                  [29531/5040, -962/35, 621/10, -4006/45, 691/8, -282/5,  2143/90,  -206/35, 363/560]]
#     assymetrical_third_derivative_parameters = [[-1,       3,      -3,         1],
#                                                 [-5/2,     9,      -12,        7,       -3/2],
#                                                 [-17/4,    71/4,   -59/2,      49/2,    -41/4,    7/4],
#                                                 [-49/8,    29,     -461/8,     62,      -307/8,   13,      -15/8],
#                                                 [-967/120, 638/15, -3929/40,   389/3,   -2545/24, 268/5,   -1849/120, 29/15],
#                                                 [-801/80,  349/6,  -18353/120, 2391/10, -1457/6,  4891/30, -561/8,    527/30, -469/240]]
#     assymetrical_fourth_derivative_parameters =[[1,       -4,       6,        -4,      1]
#                                                 [3,       -14,      26,       -24,     11,       -2]
#                                                 [35/6,    -31,      137/2,    -242/3,  107/2,    -19,      17/6]
#                                                 [28/3,    -111/2,   142,      -1219/6, 176,      -185/2,   82/3,    -7/2]
#                                                 [1069/80, -1316/15, 15289/60, -2144/5, 10993/24, -4772/15, 2803/20, -536/15, 967/240]]
#     #
#     assymetrical_derivative_parameters = [assymetrical_first_derivative_parameters,
#                                           assymetrical_second_derivative_parameters,
#                                           assymetrical_third_derivative_parameters,
#                                           assymetrical_fourth_derivative_parameters]
#     #
#     ## central difference derivative parameters
#     central_first_derivative_parameters  = [[   -1/2,      0,   1/2],
#                                             [   1/12,   -2/3,	   0,  2/3,	-1/12,],
#                                             [  -1/60,   3/20,	-3/4,	 0,	  3/4, -3/20, 1/60],
#                                             [  1/280, -4/105,	 1/5, -4/5,	    0,	 4/5, -1/5,	4/105, -1/280]]
#     central_second_derivative_parameters = [[      1,     -2,	1],				
#                                             [  -1/12,    4/3,   -5/2, 4/3,	  -1/12],			
#                                             [   1/90,  -3/20,  3/2, -49/18, 3/2,	   -3/20, 1/90],		
#                                             [ -1/560,  8/315, -1/5, 8/5,	  -205/72,	8/5,  -1/5, 8/315, -1/560]]
#     central_third_derivative_parameters  = [[   -1/2,      1,    0,       -1,    1/2],			
#                                             [    1/8,     -1,    13/8,     0,     -13/8,	1,	    -1/8],		
#                                             [ -7/240,   3/10, -169/120, 61/30, 0,	    -61/30, 169/120, -3/10, 7/240]]	
#     central_fourth_derivative_parameters = [[      1,     -4,	     6,	     -4,	      1],
# 	                                        [   -1/6,      2,	 -13/2,	   28/3,	  -13/2,	       2,	  -1/6],
# 	                                        [  7/240,   -2/5,	169/60,	-122/15,	   91/8,	 -122/15,	169/60,	-2/5,  7/240]]
#     central_fifth_derivative_parameters  = [[   -1/2,      2,	  -5/2,	      0,	    5/2,	      -2,	   1/2],
# 	                                        [    1/6,   -3/2,	  13/3,	  -29/6,	      0,	    29/6,	 -13/3,   3/2,  -1/6],
# 	                                        [-13/288,  19/36,	-87/32,	   13/2,	-323/48,	       0,	323/48, -13/2, 87/32,	-19/36,	13/288]]
#     central_sixth_derivative_parameters  = [[	   1,     -6,	    15,	    -20,	     15,	      -6,	     1],	
# 	                                        [   -1/4,      3,	   -13,	     29,	  -75/2,	      29,	   -13,	   3,  -1/4],	
# 	                                        [ 13/240, -19/24,	 87/16,	  -39/2,	  323/8,	-1023/20,	 323/8, -39/2, 87/16,	-19/24,	13/240]]
#     #
#     central_derivative_parameters = [central_first_derivative_parameters,
#                                      central_second_derivative_parameters,
#                                      central_third_derivative_parameters,
#                                      central_fourth_derivative_parameters,
#                                      central_fifth_derivative_parameters,
#                                      central_sixth_derivative_parameters]
#     #
#     function forward_difference(x, f, idx, derivative_parameters)
#         #
#         h = mean(diff(x))
#         #
#         df = 0.0
#         #
#         for i=1:lastindex(derivative_parameters)
#             df += derivative_parameters[i] * f[idx + i - 1]
#         end
#         #

#         return df / h^(d_order)
#     end
#     #  
#     function backward_difference(x, f, idx, derivative_parameters)
#         #
#         h = mean(diff(x))
#         #
#         df = 0.0
#         #
#         for i=1:lastindex(derivative_parameters)
#             df += (-1)^(d_order) * derivative_parameters[i] * f[idx - i + 1]
#         end        
#         #
#         return df / h^(d_order)
#     end
#     #   
#     function central_difference(x,f,idx,derivative_parameters)
#         #
#         h = mean(diff(x))
#         #
#         df = 0.0
#         #
#         n = Int((length(derivative_parameters)-1)/2)
#         #
#         j = 0
#         for i in -n:n
#             j += 1
#             df += derivative_parameters[j] * f[idx + i]
#         end
#         #
#         return df / h^(d_order)
#     end    
#     #
#     function max_finite_diff_orders(grid, assymetric_params, central_params)
#         n = length(grid)
#         #
#         max_forward_orders  = zeros(Int, n)
#         max_central_orders  = zeros(Int, n)
#         #
#         ## Forward and Backwards Differences
#         for i in 1:n
#             for order in 1:length(assymetric_params)
#                 if i + length(assymetric_params[order]) - 1 <= n
#                     max_forward_orders[i] = order
#                 else
#                     break
#                 end
#             end
#         end
#         #
#         ## Central differences
#         for i in 1:n
#             for order in 1:length(central_params)
#                 half_width = div(length(central_params[order]) - 1, 2)
#                 if (i - half_width >= 1) && (i + half_width <= n)
#                     max_central_orders[i] = order
#                 else
#                     break
#                 end
#             end
#         end
#         #
#         return max_forward_orders, reverse(max_forward_orders), max_central_orders
#     end
#     #
#     ## initialise derivative vector
#     dy = zeros(Float64,length(x))
#     #
#     ## maximum order for the requested derivative
#     max_len = length(assymetrical_derivative_parameters[d_order][end])
#     max_forward_order, max_backward_order, max_central_order = max_finite_diff_orders(x, assymetrical_derivative_parameters[d_order], central_derivative_parameters[d_order])
#     #
#     ## asymmetrical finite differences for the first and last 3 points up to maximum order
#     for idx=1:3
#         #
#         ## maximum order for the asymmetrical finite differences
#         F_order = max_forward_order[idx]
#         B_order = F_order #max_backward_order[idx]
#         len = Int( (length(assymetrical_derivative_parameters[d_order][F_order]) - 1))
#         #
#         ## forwards differences for first 3 points on grid
#         dy[idx] = forward_difference(x[ idx : idx + len], y[ idx : idx + len], 1, assymetrical_derivative_parameters[d_order][F_order])
#         #
#         ## backwards differences for last 3 points on grid
#         dy[end - idx + 1] = backward_difference(x[ end - idx + 1 - len : end - idx + 1 ], y[ end - idx + 1 - len : end - idx + 1 ], len + 1, assymetrical_derivative_parameters[d_order][B_order])
#     end
#     #
#     function compute_central_difference(x, y, idxs)
#         central_len = length(x[4 : end - 3])
#         #
#         if central_len % 2 == 1
#             central_idx = Int((central_len - 1)/2)
#         else
#             central_idx = Int(central_len/2)
#         end
#         #
#         dy = zeros(Float64,central_len)
#         #
#         for idx in idxs
#             C_order = Int(max_central_order[idx])
#             num_points = Int( (length(central_derivative_parameters[d_order][C_order]) - 1) / 2)
#             #
#             idx_L = 4 + idx - 1
#             idx_R = lastindex(x) - 3 - idx + 1
#             #
#             ## first half
#             dy[idx_L] = central_difference(x[idx_L - num_points: idx_L + num_points], y[idx_L - num_points: idx_L + num_points], num_points + 1, central_derivative_parameters[d_order][C_order])
#             #
#             ## last half
#             dy[idx_R] = central_difference(x[idx_R - num_points: idx_R + num_points], y[idx_R - num_points: idx_R + num_points], num_points + 1, central_derivative_parameters[d_order][C_order])
#         end
#         #
#         if central_len % 2 == 1
#             idx_C      = start_idx + 3 + central_idx
#             C_order    = Int(max_central_order[idx_C])
#             num_points = Int( (length(central_derivative_parameters[d_order][C_order]) - 1) / 2)
#             dy[idx_C]  = central_difference(x[idx_C - num_points: idx_C +  num_points], y[idx_C -  num_points: idx_C +  num_points],  num_points + 1,  central_derivative_parameters[d_order][C_order])
#         end
#         #
#         return dy
#     end
#     #
#     ## multi-threading
#     # if Threads.nthreads() > 1
#     #     idxs = collect(4:lastindex(x)-3)
#     #     num_threads = min(Threads.nthreads(), length(idxs))
#     #     chunk_size = ceil(Int, length(idxs) / num_threads)
#     #     #
#         # @threads for i in 1:num_threads
#         #     start_idx = idxs[(i - 1) * chunk_size + 1]
#         #     end_idx = Int(min(i * chunk_size, length(idxs)))
#         #     end_idx = idxs[end_idx]
#         #     #
#         #     if start_idx <= end_idx
#         #         for idx=start_idx:end_idx
#         #             local C_order, num_points
#         #             #
#         #             C_order = Int(max_central_order[idx])
#         #             num_points = Int( (length(central_derivative_parameters[d_order][C_order]) - 1) / 2)
#         #             #
#         #             dy[idx] = central_difference(x[idx - num_points: idx + num_points], y[idx - num_points: idx + num_points], num_points + 1, central_derivative_parameters[d_order][C_order])
#         #         end
#         #     end
#         # end
#     # else
#     #
#     ## central differences for middle points
#     central_len = length(x[4:lastindex(x)-3])
#     if central_len % 2 == 1
#         central_idx = Int((central_len - 1)/2)
#     else
#         central_idx = Int(central_len/2)
#     end
#     #
#     for idx=1:central_idx
#         C_order = Int(max_central_order[4 + idx - 1])
#         num_points = Int( (length(central_derivative_parameters[d_order][C_order]) - 1) / 2)
#         #
#         idx_L = 4 + idx - 1
#         idx_R = lastindex(x) - 3 - idx + 1
#         #
#         ## first half
#         dy[idx_L] = central_difference(x[idx_L - num_points: idx_L + num_points], y[idx_L - num_points: idx_L + num_points], num_points + 1, central_derivative_parameters[d_order][C_order])
#         #
#         ## last half
#         dy[idx_R] = central_difference(x[idx_R - num_points: idx_R + num_points], y[idx_R - num_points: idx_R + num_points], num_points + 1, central_derivative_parameters[d_order][C_order])
#     end

#     if central_len % 2 == 1
#         idx_C      = 4 + central_idx
#         C_order    = Int(max_central_order[idx_C])
#         num_points = Int( (length(central_derivative_parameters[d_order][C_order]) - 1) / 2)
#         dy[idx_C]  = central_difference(x[idx_C - num_points: idx_C +  num_points], y[idx_C -  num_points: idx_C +  num_points],  num_points + 1,  central_derivative_parameters[d_order][C_order])
#     end
#     # end
#     #
#     return dy
# end
# #
# function adiabatic_HO(r, V1d, V2d, gamma, re)
#     #
#     ##
#     beta = pi/4 .+ 0.5 .* atan.((r .- re) ./ gamma)
#     W12 = 0.5 .* gamma .* ( (r.-re).^2 .+ gamma^2).^(-1)
#     #
#     ∂Vd = V2d .- V1d
#     #
#     DC = 0.5 .* tan.(2 .* beta) .* ∂Vd
#     #
#     ##
#     V1a = (V1d .+ V2d) .- sqrt.(∂Vd.^2 + 4 .* DC.^2) #V1d .* cos.(beta).^2 .+ V2d .* sin.(beta).^2
#     V2a = (V1d .+ V2d) .+ sqrt.(∂Vd.^2 + 4 .* DC.^2) #V2d .* cos.(beta).^2 .+ V1d .* sin.(beta).^2
#     #
#     return [V1a, V2a], W12, beta
# end
# #
# function repulsive(r,C3,C6)
#     return C3*r^(-3) + C6*r^(-6)
# end



# y = state_v_wav[1,1,:]

# d2y = FiniteDifference(collect(r),y,2)

# println("Ke_11 from sinc = ",T[1,1]," from integral = ",simps(y .* d2y, r[1],r[end]))


# dim = Ne*vmax

# T_v = zeros(Float64, dim, dim)

# for state_i = 1 : Ne
#     #
#     for v_i=1:vmax
#         for v_j=v_i:vmax
#             #
#             ψ_statei_vi = state_v_wav[state_i,v_i,:]
#             ψ_statej_vj = state_v_wav[state_j,v_j,:]
#             #
#             T_v = simps(state_v_wav[state_i,v_i,:] .* FiniteDifference(r,state_v_wav[state_j,v_j,:],2), r[1],r[end])

#     for  state_j = state_i + 1 : Ne















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
