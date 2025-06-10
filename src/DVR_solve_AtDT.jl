using LegendrePolynomials
using PyPlot
using FastGaussQuadrature
using LinearAlgebra
using Dierckx
using Statistics
#
function cumtrapz(X::T, Y::T) where {T <: AbstractVector}
    # Check matching vector length
    @assert length(X) == length(Y)
    # Initialize Output
    out = similar(X)
    out[1] = 0
    # Iterate over arrays
    for i in 2:length(X)
      out[i] = out[i-1] + 0.5*(X[i] - X[i-1])*(Y[i] + Y[i-1])
    end
    # Return output
    out
end
#
function plot_heatmap_pyplot_julia_highres(matrix)
    """
    Plots a matrix as a heatmap with high color resolution.

    Args:
        matrix: The matrix to plot.
        title_str: (Optional) The title of the plot. Defaults to "High-Resolution Heatmap".
    """

    vmin = minimum(matrix)
    vmax = maximum(matrix)
    plt.figure()
    plt.imshow(matrix, cmap="plasma", origin="lower", aspect="auto", vmin=vmin, vmax=vmax)
    colorbar()
    xlabel("Column Index")
    ylabel("Row Index")
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
function derivative_sinc(x::Float64, xj::Float64, h::Float64)::Float64
    # Avoid division by zero: handle x == x_j separately
    if abs(x - xj) < 1e-10
        return 0.0
    else
        return ( pi^2 * (x-xj) * cos(pi * (x - xj) / h) - h*pi*sin(pi * (x - xj) / h) ) / (pi * (x - xj))^2
    end
end
#
@views function simps(f::Function, x::AbstractRange) #A range with odd number
    h = step(x)
    I= h/3*(f(x[1])+2*sum(f,x[3:2:end-2])+4*sum(f,x[2:2:end-1])+f(x[end]))
    return I
end
#
# from the function and range definition
function simps(f::Function, a::Real, b::Real, n::Integer) #n as an even number
    return simps(f, range(a, b, length=n+1))
end
#
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
function lorentzian_reNorm_N2(r, N, w, r0)      
    return N*0.5*(w/((r - r0)^2 + w^2))
end
#

Ngrid = 10001
r = LinRange(0,10,Ngrid)
h = r[2]-r[1]



p12 = [1,0.02699694, 2.31901939]
p13 = [0.095, 0.14, 2.34]
p23 = [-1,0.0477113 ,2.10810043] 

W12 = lorentzian_reNorm_N2.(r,p12...)
W13 = lorentzian_reNorm_N2.(r,p13...)
W23 = lorentzian_reNorm_N2.(r,p23...)
#
W_ = map(x -> [ 0 W12[x] W13[x] ; -W12[x] 0 W23[x] ; -W13[x] -W23[x] 0], collect(1:size(r)[1]))
#
W = zeros(Float64,Ngrid,dim,dim)
for i=1:dim
    for j=i+1:dim
        W[:,i,j] .=  [w[i,j] for w in W_]
        W[:,j,i] .= -[w[i,j] for w in W_]
    end
end

norm = sqrt.(W12.^2 .+ W13.^2 .+ W23.^2)

cdf = cumtrapz(collect(r),norm)


# plt.figure()
# plt.plot(cdf,W12./norm)
# plt.plot(cdf,W13./norm)
# plt.plot(cdf,W23./norm)

cdf_spline = Spline1D(r,cdf)

# cdf = cdf ./cdf[end]        # normalise to 1
#
## find the inverse cdf and spline it
cdf_inverted_spline = Spline1D(cdf,r)
Ngrid = 1001
#
## now compute the inverse transform sample of the cdf to yield non-unifrom nuclear configuration grid
r_grid = cdf_inverted_spline(LinRange(0,cdf[end],Ngrid))
#
## re compute NACs
W12 = lorentzian_reNorm_N2.(r_grid,p12...)
W13 = lorentzian_reNorm_N2.(r_grid,p13...)
W23 = lorentzian_reNorm_N2.(r_grid,p23...)
#
norm = sqrt.(W12.^2 .+ W13.^2 .+ W23.^2)
#
W_ = map(x -> [ 0 W12[x] W13[x] ; -W12[x] 0 W23[x] ; -W13[x] -W23[x] 0], collect(1:size(r_grid)[1]))
#
W = zeros(Float64,length(r_grid),dim,dim)
for i=1:dim
    for j=i+1:dim
        W[:,i,j] .=  [w[i,j] for w in W_] ./ norm
        W[:,j,i] .= -[w[i,j] for w in W_] ./ norm
    end
end


plt.figure()
plt.plot(r_grid,W[:,1,2])
plt.plot(r_grid,W[:,1,3])
plt.plot(r_grid,W[:,2,3])



X = collect(LinRange(0,cdf[end],Ngrid))
h = X[2]-X[1]



# y = sinc_function.(r,r[500],diff(r)[1])

# dy = derivative_sinc.(r,r[500],diff(r)[1])

# plt.figure()
# plt.plot(r,y)
# plt.yscale("log")

# plt.figure()
# plt.plot(r,dy)

dim = 3
N = dim*Ngrid
M = zeros(Float64, Ngrid, Ngrid)
#
for i=1:Ngrid
    for j=1:Ngrid
        if i==j
            M[i,j] = 0.0
        else
            M[i,j] = (-1)^(i-j)/((i-j))
        end
    end
end
#
I = zeros(Float64,dim,dim)
for d=1:dim
    # p = zeros(Float64,dim,dim)
    # for i=1:dim
    #     p[i,d] = 1.0
    # end
    # #
    # push!(P,kron(p,M))
    I[d,d] = 1.0
end

#
Y = kron(I,M)
println(1)
plot_heatmap_pyplot_julia_highres(Y)


Wdvr = Matrix{Float64}[]
B = zeros(Float64,N,N)

for b=1:dim
    for g=1:dim
        if b != g
            Wbg = W[:,b,g]
            #
            for idx=1:Ngrid
                i = Int((b-1)*Ngrid + idx)
                f = Int((g-1)*Ngrid + idx)
                #
                B[i,f] = Wbg[idx]*h
            end
        end
    end
end
println(2)
# plot_heatmap_pyplot_julia_highres(B)
#
S = Y + B
#
c = nullspace(S)
#
U_vec = []

for d=1:dim
    push!(U_vec,c[(d-1)*Ngrid+1:d*Ngrid])
    # println((d-1)*Ngrid+1," ",d*Ngrid)
end

plt.figure()
plt.plot(X,U_vec[1])
plt.plot(X,U_vec[2])
plt.plot(X,U_vec[3])







# # for g=1:dim
# #     B_ = zeros(Float64,Ngrid,Ngrid)
# #     if g!=b
# #         Wbg = W[:,b,g]
# #         #
# #         for i=1:Ngrid
# #             B_[i,i] = Wbg[i]
# #         end
# #     end
# #     #
# #     push!(Wdvr,B_)
# # end
# # #
# # B = BlockDiagonal(Wdvr)

# # x = zeros(Float64,Ngrid,Ngrid)
# # for i=1:Ngrid
# #     for j=1:Ngrid
# #         x[i,j] = h*sinc_function(r[i],r[i],h)*derivative_sinc(r[i],r[j],h)
# #     end
# # end

# # z= zeros(Float64,Ngrid,Ngrid)


# # for i=1:Ngrid
# #     for j=1:Ngrid
# #         z[i,j] = h * sum(sinc_function.(r,r[i],h) .* sinc_function.(r,r[j],h) .* W12)
# #         println(sinc_function.(r,r[1],h)," ",sinc_function.(r,r[j],h)," ",W12,sum(sinc_function.(r,r[i],h) .* sinc_function.(r,r[j],h) .* W12))
# #     end
# # end









# # for i=1:Ngrid
# #     for j=1:Ngrid
# #         if i==j
# #             M[i,j] = 0.0
# #         else
# #             M[i,j] = (-1)^(i-j)/(i-j)
# #         end
# #     end
# # end

