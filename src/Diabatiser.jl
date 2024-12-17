########################################################################################################################
# 
# ABOUT: THIS MODULE CONTAINS ALL OF THE FUNCTIONS REQUIRED TO DIABATISE A SET(S) OF ADIABATICALLY INTERACTING STATES  #
#        AND RETURN THE CORESPONDING SMOOTH DIABATIC POTENTIALS (PECs). THIS MODULE USES THE PROPERTY OF THE           #
#        ELECTRONIC WAVEFUNCTIONS SHOULD BE SQUARE INTEGRABLE AND CONTINUOUS SUCH THAT THEIR PECs ARE ALSO SMOOTH IN   #
#        THE DIABATIC REPRESENTATION, SO NON-ADIABATIC-COUPLINGS (NACs) ARE FITTED ON THE BASIS THAT THE DIABATISED    #
#        PECs SHOULD BE SMOOTH.                                                                                        #
#           NOT ONLY THE 2X2 PROBLEM CAN BE SOLVED, BUT THE ROUTINES HERE CAN DIABATISE AN ARBITRARY DIMENSIONED       #
#        CROSSING SYSTEM.                                                                                              #
#                                                                                                                      #
# FOR MORE DETAILS ON THE METHODOLGY PLEASE REFER TO DOI: 10.1039/D2CP03051A OR CONTACT RYAN BRADY AT:                 #
#                                                                                                                      #
#                                                 ryan.brady.17@ucl.ac.uk                                              #
# LIBRARY IMPORT:                                                                                                      #
# using LinearAlgebra
# using Optim
#                                                                                                                      #
# SETTING ENVIRONMENT VARIABLES                                                                                        #
# ENV["JULIA_NUM_THREADS"] = 7
#                                                                                                                      #
# EXPORTING FUNCTIONS                                                                                                  #
# export lorentzian, mixAng_lo, mixAng_avg, mixAng_la, laplace, diabatise, adiabatise, fit_diabat, diabatise_ao
#                                                                                                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Exceptions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
## No found crossings
# mutable struct NoCrossing <: Exception
#         var::Vector{Number}
# end
# NoCrossingMessage = "No Crossings where found. Try decreasing the 'peak_prominence' parameter of 'Diabatise_pipe'. 
#                   If this doesn't work, please try changing the 'min_peak_distance' and 'thresh' parameters. 
#                   See: https://juliapackages.com/p/findpeaks for more details. \n
#                   Peak Max, Min, average are: "

# Base.showerror(io::IO, e::NoCrossing) = print(io, NoCrossingMessage, e.var, "\n")
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Functional forms for NACs and Mixing Angles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
function display_progress(i,j,bar_length,current,total)
    text = "Diabatising Pair "*string(i)*","*string(j)
    #
    percent_complete = round(current / total * 100)
    #
    num_bar_symbols = Int(round(current / total * bar_length))  # Convert to integer
    bar = "" * repeat("█" , num_bar_symbols) * repeat(" ", bar_length - num_bar_symbols) * ""
    
    # Use carriage return to overwrite the previous line
    print("\r"*text*": $percent_complete% $bar ($current/$total)")
    flush(stdout)
end
#
"""
ABOUT:  get_loss() takes pairs of adiabatically interacting states with their crossing (NAC) parameters and will perform
        the diabatisation in the 2-state approximation and return the smoothness of the diabatic property (frobenius
        norm of second derivative matrix summed).
------------------------------------------------------------------------------------------------------------------------
INPUT:
        1) r: bond length,
        2) a: diagonal 2x2 matrix of adiabatically interacting states,
        3) f: functional form for the NAC,
        4) p: function f(p) parameters.
------------------------------------------------------------------------------------------------------------------------
OUTPUT:
        1) loss: the sum of the absolute second order derivatives of the upper diabatic state.
"""
function get_loss(r, Lval, units, sub_type, factor, adiabats, f, fixParams, p_init, p_bounds, p)
        #
        ## finite difference (high order accuracy) for matrices
        function FiniteDiff_MatDerivative(x, M, dim, d_order)
            dM = zeros(Float64,length(x),dim,dim)
            #
            ## extract elements of U, spline them, form splined matrix and derivative
            for i=1:dim
                for j=1:dim
                    #
                    Mij = [u[i,j] for u in M]
                    #
                    dMij = FiniteDifference(x,Mij,d_order)
                    #
                    dM[:,i,j] = dMij
                end
            end
            #
            return [dM[idx,:,:] for idx=1:lastindex(x)]
        end
        #
        ## smoothness of the dynamics induced into diabats by U (omitting adiabats)
        function SmoothnessByU(r,U,Va,dim)
            Vd = adjoint.(U) .* Va .* U
            #
            d2Vd = Vector{SMatrix{dim, dim, Float64}}(FiniteDiff_MatDerivative(r,Vd,dim,2))
            #
            SmoothSum = zeros(Float64,dim)
            #
            for x=1:lastindex(r)
                #
                # M = d2U[x]'*Va[x]*U[x] + 2*dU[x]'*Va[x]*dU[x] + U[x]'*Va[x]*d2U[x]
                M = d2Vd[x]
                #
                for i=1:dim
                    for j=i:dim
                        SmoothSum[i] += M[i,j]^2
                    end
                end
            end
        #     println(sum(SmoothSum))
            #
            return sqrt(sum(SmoothSum))
        end
        #
        function SmoothnessByU_old(r,U,Va,dim)
                #
                dU  = Vector{SMatrix{dim, dim, Float64}}(FiniteDiff_MatDerivative(r,U,dim,1))
                d2U = Vector{SMatrix{dim, dim, Float64}}(FiniteDiff_MatDerivative(r,U,dim,2))
                #
                U  = Vector{SMatrix{dim, dim, Float64}}(U)
                Va = Vector{SMatrix{dim, dim, Float64}}(Va)
                #
                SmoothSum = zeros(Float64,dim)
                #
                for x=1:lastindex(r)
                    #
                    M = d2U[x]'*Va[x]*U[x] + 2*dU[x]'*Va[x]*dU[x] + U[x]'*Va[x]*d2U[x]
                    #
                    for i=1:dim
                        for j=i:dim
                            SmoothSum[i] += M[i,j]^2
                        end
                    end
                end
                #
                return sqrt(sum(SmoothSum))
        end
        #
        for (i,par) in enumerate(p)
                if fixParams[i] == 0.0
                        p[i] = p_init[i]
                else
                        p[i] = par #round(par, sigdigits=6)
                        #
                        if (p[i]<p_bounds[i][1])|(p[i]>p_bounds[i][2])
                                return 1e100
                        end
                end
        end
        #
        ## compute mixing angle
        if in(f, ["lorentzian","laplacian","gaussian","lor_lap"])
                #
                mixing_angle_func = eval(Symbol("mixAng_"*f))
                #
                mixing_angle = mixing_angle_func.(r,p...) .* factor
        else
                X, NACij = ComputeProperty_viaParameters(r, f, Lval, p, "NAC", units, sub_type, factor)
                mixing_angle = map(i->trapz(r[1:i],NACij[1:i]), collect(1:size(r)[1]))
        end
        #
        ## compute the 2x2 diabatising Matrix
        U = map(beta -> [cos(beta) -sin(beta) ; sin(beta) cos(beta)], mixing_angle)
        #
        ## compute smoothness
        cost = SmoothnessByU(r,U,adiabats,2)
        #
        return round(cost, sigdigits=6)
end
#
function check_symmetry_and_states(i,j)
        symmetry_conditions = [Potential[i].mult     == Potential[j].mult,
                               Potential[i].lambda   == Potential[j].lambda,
                               Potential[i].symmetry == Potential[j].symmetry,
                               i in Calculation["method"].states,
                               j in Calculation["method"].states]
        #
        return all(symmetry_conditions)
end
#
function compute_mixing_angle(r, f, key, p)
        if in(f, ["lorentzian","laplacian","gaussian","lor_lap"])
                #
                mixing_angle_func = eval(Symbol("mixAng_"*f))
                #
                mixing_angle = mixing_angle_func.(r,p...) .* NonAdiabaticCoupling[key].factor
        else
                X, NACij = ComputeProperty_viaParameters(r, f, NonAdiabaticCoupling[key].Lval, p, "NAC", NonAdiabaticCoupling[key].units, NonAdiabaticCoupling[key].sub_type,  NonAdiabaticCoupling[key].factor)
                #
                mixing_angle = map(i->trapz(X[1:i],NACij[1:i]), collect(1:size(r)[1]))
        end
        #
        return mixing_angle
end
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DIABATISING ROUTINES  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
## 2-state approximation
# function fit_multi_diabat_2stateApprox(r, a::Array; precision = 1e-6)
#         options = Optim.Options(show_trace = true)
#         #
#         D  = Int(size(a)[2])
#         #
#         ## initialise the adiabatic to diabatic transformation
#         U = zeros(Float64, size(a)[1], D, D)
#         #
#         for i=1:D
#                 U[:,i,i] .= 1.0
#         end
#         #
#         ## loop over different NACs
#         xidx = 0
#         for key in keys(NonAdiabaticCoupling)
#                 xidx +=1
#                 #
#                 ## bra and ket labels
#                 i = floor(Int64, key[1])
#                 j = floor(Int64, key[2])
#                 #
#                 display_progress(i,j, 16, xidx, length(keys(NonAdiabaticCoupling)))
#                 print("\n")
#                 #
#                 ## extract fitting flags, i.e. turn of parameter variation in fit
#                 p_excludeFromFit = NonAdiabaticCoupling[key].fit
#                 #
#                 ## extract guesses for parameters
#                 p_guess = NonAdiabaticCoupling[key].Rval
#                 #
#                 ## determine the functional form to fit
#                 func = NonAdiabaticCoupling[key].type
#                 #
#                 ## Generate 2x2 problem
#                 a_tmp = Array{Float64}(undef, size(a)[1],2,2) 
#                 for k=1:size(a)[1]
#                         a_tmp[k,:,:] = Diagonal(Array([a[k,i,i],a[k,j,j]])) 
#                 end
#                 a_tmp = [a_tmp[idx,:,:] for idx=1:lastindex(r)]
#                 #
#                 ## optimise NAC parameters to create the smoothest PECs for states i,j
#                 f = func
#                 p =  p_guess
#                 #
#                 ##
#                 mixing_angle_ij = compute_mixing_angle(r, f, key, p)
#                 #
#                 c = cos.(mixing_angle_ij)
#                 s = sin.(mixing_angle_ij)
#                 #
#                 U[:,i,i] .= c
#                 U[:,j,j] .= c
#                 U[:,i,j] .= -s
#                 U[:,j,i] .= s
#                 #
#                 ## if no parameters are to be fit then skip fitting step
#                 if (any(p_excludeFromFit .== 1))&(check_symmetry_and_states(i,j))
#                         #
#                         ## extract parameter bounds
#                         p_bounds = NonAdiabaticCoupling[key].bounds
#                         #
#                         ## obtain parameters for the ComputePropert_viaParameters function
#                         Lval = NonAdiabaticCoupling[key].Lval
#                         units = NonAdiabaticCoupling[key].units
#                         sub_type = NonAdiabaticCoupling[key].sub_type
#                         factor = NonAdiabaticCoupling[key].factor
#                         #
#                         ## perform optimization
#                         o_ = optimize(p -> get_loss(r, Lval, units, sub_type, factor, a_tmp, func, p_excludeFromFit, p_guess, p_bounds, p), [p_guess...],options) 
#                         optimisedParameters = Optim.minimizer(o_)
#                         #
#                         NonAdiabaticCoupling[key].fitted_parameters .= optimisedParameters
#                         #
#                         ## compute mixing angles
#                         mixing_angle_ij = compute_mixing_angle(r,f,key,optimisedParameters)
#                         #
#                         c = cos.(mixing_angle)
#                         s = sin.(mixing_angle)
#                         #
#                         U[:,i,i] .= c
#                         U[:,j,j] .= c
#                         U[:,i,j] .= -s
#                         U[:,j,i] .= s
#                 end
#         end
#         #
#         return U
# end
function fit_multi_diabat_2stateApprox(r, a::Array; precision = 1e-6)
    options = Optim.Options(show_trace = true)
    #
    D  = Int(size(a)[2])
    #
    ## initialise the adiabatic to diabatic transformation
    U = zeros(Float64, size(a)[1], D, D)
    #
    for i=1:D
            U[:,i,i] .= 1.0
    end
    #
    ## loop over different NACs
    xidx = 0
    for key in keys(NonAdiabaticCoupling)
            xidx +=1
            #
            ## bra and ket labels
            i = floor(Int64, key[1])
            j = floor(Int64, key[2])
            #
            if (i in Calculation["method"].states)&(j in Calculation["method"].states)
                #
                display_progress(i,j, 16, xidx, length(keys(NonAdiabaticCoupling)))
                print("\n")
                #
                ## extract fitting flags, i.e. turn of parameter variation in fit
                p_excludeFromFit = NonAdiabaticCoupling[key].fit
                #
                ## extract guesses for parameters
                p_guess = NonAdiabaticCoupling[key].Rval
                #
                ## determine the functional form to fit
                func = NonAdiabaticCoupling[key].type
                #
                ## Generate 2x2 problem
                a_tmp = Array{Float64}(undef, size(a)[1],2,2) 
                for k=1:size(a)[1]
                        a_tmp[k,:,:] = Diagonal(Array([a[k,i,i],a[k,j,j]])) 
                end
                a_tmp = [a_tmp[idx,:,:] for idx=1:lastindex(r)]
                #
                ## optimise NAC parameters to create the smoothest PECs for states i,j
                f = func
                p =  p_guess
                #
                ##
                mixing_angle_ij = compute_mixing_angle(r, f, key, p)
                #
                c = cos.(mixing_angle_ij)
                s = sin.(mixing_angle_ij)
                #
                U[:,i,i] .= c
                U[:,j,j] .= c
                U[:,i,j] .= -s
                U[:,j,i] .= s
                #
                ## if no parameters are to be fit then skip fitting step
                if (any(p_excludeFromFit .== 1))&(check_symmetry_and_states(i,j))
                        #
                        ## extract parameter bounds
                        p_bounds = NonAdiabaticCoupling[key].bounds
                        #
                        ## obtain parameters for the ComputePropert_viaParameters function
                        Lval = NonAdiabaticCoupling[key].Lval
                        units = NonAdiabaticCoupling[key].units
                        sub_type = NonAdiabaticCoupling[key].sub_type
                        factor = NonAdiabaticCoupling[key].factor
                        #
                        ## perform optimization
                        o_ = optimize(p -> get_loss(r, Lval, units, sub_type, factor, a_tmp, func, p_excludeFromFit, p_guess, p_bounds, p), [p_guess...],options) 
                        optimisedParameters = Optim.minimizer(o_)
                        #
                        NonAdiabaticCoupling[key].fitted_parameters .= optimisedParameters
                        #
                        ## compute mixing angles
                        mixing_angle_ij = compute_mixing_angle(r,f,key,optimisedParameters)
                        #
                        c = cos.(mixing_angle_ij)
                        s = sin.(mixing_angle_ij)
                        #
                        U[:,i,i] .= c
                        U[:,j,j] .= c
                        U[:,i,j] .= -s
                        U[:,j,i] .= s
                end
            end
    end
    #
    return U
end
#
## N-state
#
function dynamic_nac(x)
    #
    ## initialse the evolving NAC matrix
    Wi = zeros(Float64,dim,dim)
    #
    for key in keys(NonAdiabaticCoupling)
        local _, NAC
        _, NAC = ComputeProperty(NonAdiabaticCoupling[key], custom_grid = [x], evolution_grid = false)
        Wi[key[1],key[2]] =   NAC[1] .* NonAdiabaticCoupling[key].factor
        Wi[key[2],key[1]] = - NAC[1] .* NonAdiabaticCoupling[key].factor
    end
    #
    return Wi
end
#
function adaptive_grid(error, x; µ = 1e-12)
    if (x[end] < 0)|(x[end]>10)
        ∂x_tol = 1
    else
        rang = Calculation["grid"].range
        ∂x_tol = 1e-3 #(rang[2]-rang[1])/Calculation["grid"].npoints
    end
    # println(error[end-1:end],x[end-1:end])
    #
    ## push the grid point in opposite direction of the derivative of error
    dE = (error[end] - error[end-1])/((x[end] - x[end - 1]))
    #
    ## compute approximate new position
    ∂x = min((µ - error[end])/(dE), ∂x_tol)
    x_current = x[end] + ∂x
    return x_current
end
#
function KE_error(U,dU,W; tol = 1e-12)
    K = dU'*dU 
    #
    ## residual KE
    Ke = K - U'*W*W*U + dU'*W*U - U'*W*dU
    #
    return min(frobeniusNorm(Ke*2), tol) #Calculation["method"].KE_factor
end
#
function Kangle(K)
        alpha = []
        for i=1:size(K)[1]
            for j=i+1:size(K)[2]
                push!(alpha,K[i,j]^2)
            end
        end
        return sqrt(sum(alpha))
end
#
function EvolutionOperator(rf,ri,Wf,Wi,dim,Identity)
    #
    ## form the matrix of CDFs (betas) by compute of the integral matrix, approx trapzoidal rule
    K = 0.5*(rf-ri)*(Wf+Wi)
    #
    ## find the exponential of these matrices
    if dim == 3                                            # use rodriguez formular
        alpha = Kangle(K)                                  # find the corresponding combined mixing half angle / cdf
        #
        ## compute evolution matrix by rodriguez rotation formular
        evoMat = Identity - (sin(alpha)/alpha)*K + ((1-cos(alpha))/alpha^2)*K^2
    elseif dim > 3                                         # use diagonalisation/Schur decomposition
        evoMat = exp(-K)
    end
    #
    return evoMat, K
end
#
function Forward_Evolution(rgrid, W, dim, BoundaryCondition)
    #
    ## compute the identity matrix
    Identity = zeros(Float64,dim,dim)
    for i=1:dim
        Identity[i,i] = 1.0
    end
    #
    ## number of points
    N = length(rgrid)
    #
    ## now evolve the matrix from unit matrix
    U = zeros(Float64,N,dim,dim)
    U[1,:,:] = BoundaryCondition
    #
    for idx=2:lastindex(rgrid)
        #
        ## compute the initial evolution prior to regularisation
        E, W_int = EvolutionOperator(rgrid[idx],rgrid[idx-1],W[idx,:,:],W[idx-1,:,:],dim,Identity)
        #
        ## perform local evolution
        U[idx,:,:] = E*U[idx-1,:,:]
    end
    #
    ## compute new NACs
    dU = FiniteDiff_MatDerivative(rgrid,U,dim, 1) # computes derivative of matrix
    UdU = map(i -> U[i,:,:]*dU[i,:,:]',collect(1:size(rgrid)[1]))
    #
    return U, dU, UdU
end
#
function Backward_Evolution(rgrid, W, dim, BoundaryCondition)
    #
    ## compute the identity matrix
    Identity = zeros(Float64,dim,dim)
    for i=1:dim
        Identity[i,i] = 1.0
    end
    #
    ## number of points
    N = length(rgrid)
    #
    ## now evolve the matrix from unit matrix
    U = zeros(Float64,N,dim,dim)
    U[end,:,:] = BoundaryCondition
    #
    for idx in reverse(1:N-1)
        #
        ## compute the initial evolution prior to regularisation
        E, W_int = EvolutionOperator(rgrid[idx],rgrid[idx+1],W[idx,:,:],W[idx+1,:,:],dim,Identity)
        #
        ## perform local evolution
        U[idx,:,:] = E*U[idx+1,:,:]
    end
    #
    ## compute new NACs
    dU = FiniteDiff_MatDerivative(rgrid,U,dim, 1) # computes derivative of matrix
    UdU = map(i -> U[i,:,:]*dU[i,:,:]',collect(1:size(rgrid)[1]))
    #
    return U, dU, UdU
end
#
function find_closest_permutation(A,dim)
    Q = zeros(Float64,dim,dim)
    #
    for i=1:dim
        j = findall(abs.(A[:,i]) .== maximum(abs.(A[:,i])))[1]
        #
        if A[j,i] < 0 
            sign = -1.0
        else
            sign = 1.0
        end
        #
        Q[j,i] = sign
    end
    #
    return Q
end
#
function regularise_old(rgrid, Uf, Ub, Va, dim, order, region; DENSE_GRID_OBJECTS = false, µ = 0.1)
    rgrid = collect(rgrid)
    #                  
    ######################### FUNCTIONS FOR SUBROUTINE #########################
    #
    ## switching generator matrices with #dim sigmoids
    function matMix_const_gamma(x, A, B, r0_vec, gamma_vec)
        #
        ## initialise matrix C
        C = zeros(Float64,dim,dim)
        #
        ## assume row stacking indices
        k = 0
        #
        for i=1:dim
            for j=i+1:dim
                #
                k += 1
                #
                f = sigmoid(x, gamma_vec[k], r0_vec[k])
                #
                C[i,j] = f * B[i,j] + (1 - f) * A[i,j] 
                C[j,i] = f * B[j,i] + (1 - f) * A[j,i]
            end
        end
        #
        return C
    end
    #
    ## cost function for switching functions with constant slope (gamma) params
    function cost_init_gamma(rgrid, Va, Kf, Kb, dim, r0, gamma, gamma_bounds, r0_bounds)
        #
        ## compute adaptive gamma bounds
        # gamma_bounds = []
        # k = 0
        # for i=1:dim
        #     for j=i+1:dim
        #         k += 1
        #         g_thresh_left  = gamma_sigmoid(region[1],   µ, r0[k])
        #         g_thresh_right = gamma_sigmoid(region[2], 1-µ, r0[k])
        #         #
        #         push!(gamma_bounds, minimum([g_thresh_left, g_thresh_right]))
        #         println([g_thresh_left, g_thresh_right], r0[k])
        #     end
        # end
        # # = []
        r0 = round.(r0, sigdigits=6) 
        r0_bool = []
        k = 0
        for i=1:dim
            for j=i+1:dim
                k+=1
                #
                push!(r0_bool, r0_bounds[k][1] <= r0[k] <= r0_bounds[k][2])
            end
        end
        # 
        ∂g = gamma .- gamma_bounds
        # println(gamma_bounds,"  ",r0)
        #
        gamma = round.(gamma, sigdigits=2) 
        #
        if any(∂g .< 0.0)
            return 1e100
        elseif !all(r0_bool) #any([a < rgrid[1] || a > rgrid[end] for a in r0])
            return 1e100
        end
        #
        ## compute mixing of Kf and Kb
        K  = map( x -> matMix_const_gamma(rgrid[x], Kf[x], Kb[x], r0, gamma), collect(1:lastindex(rgrid)))
        #
        ## compute the corresponding unitary matrix
        U  = real(exp.(K))
        #
        ## compute smoothness of the diabatic curves induced by U (ignore adiabatic smoothness/derivatives)
        smoothness =  SmoothnessByU(rgrid,U,Va,dim)
        #
        return smoothness
    end
    #
    ## sigmoid with r-dependent gamma
    function sigmoid_variable_gamma(r,gamma,r0)
        return map(x -> 1/(1+exp(-gamma[x]*(r[x]-r0))), collect(1:lastindex(r)))
    end
    #
    ## morphing function
    function polynom_decay_24(r, re, B; Binf = 1, beta2 = 0.1, beta4 = 0.02, p = 4)
        #
        ## Surkus variable
        y = (r^p-re^p)/(r^p+re^p)
        #
        ## reduced variable
        z = (r-re) * exp(-beta2*(r-re)^2-beta4*(r-re)^4) # sigmoid_square(r, region[1], region[2], re)
        #
        expansion = sum([B[k] * z^(k-1) * (1-y) for k=1:size(B)[1]]) + Binf * y
        #
        return expansion
    end
    #
    ## morph constant to a curve via multiplication with f =  polynom_decay_24
    function morph_const_pd24(r, fconst, re, B; Binf = 1, beta2 = 0.1, beta4 = 0.02, p = 4)
        #
        return fconst .* map(x -> polynom_decay_24(x, re, B, Binf = Binf, beta2 = beta2, beta4 = beta4, p = p), r)
    end
    #
    ## finite difference (high order accuracy) for matrices
    function FiniteDiff_MatDerivative(x, M, dim, d_order)
        dM = zeros(Float64,length(x),dim,dim)
        #
        ## extract elements of U, spline them, form splined matrix and derivative
        for i=1:dim
            for j=1:dim
                #
                Mij = [u[i,j] for u in M]
                #
                dMij = FiniteDifference(x,Mij,d_order)
                #
                dM[:,i,j] = dMij
            end
        end
        #
        return [dM[idx,:,:] for idx=1:lastindex(x)]
    end
    #
    ## smoothness of the dynamics induced into diabats by U (omitting adiabats)
    function SmoothnessByU(r,U,Va,dim)
        #
        # if Threads.nthreads() >=2
        #     #
        #     derivs = Vector{Vector{SMatrix{dim, dim, Float64}}}(undef, 2)
        #     #
        #     @threads for i in 1:2
        #         derivs[i] = Vector{SMatrix{dim, dim, Float64}}(FiniteDiff_MatDerivative(r,U,dim,i))
        #     end
        #     #
        #     dU, d2U = derivs[1], derivs[2]
        # else
        dU  = Vector{SMatrix{dim, dim, Float64}}(FiniteDiff_MatDerivative(r,U,dim,1))
        d2U = Vector{SMatrix{dim, dim, Float64}}(FiniteDiff_MatDerivative(r,U,dim,2))
        # end
        #
        U  = Vector{SMatrix{dim, dim, Float64}}(U)
        Va = Vector{SMatrix{dim, dim, Float64}}(Va)
        #
        # n = length(r)
        # num_threads = min(Threads.nthreads(), n)
        # chunk_size = ceil(Int, n / num_threads)
        
        # SmoothSum = zeros(Float64,dim)
    
        # @threads for i in 1:num_threads
        #     start_idx = (i - 1) * chunk_size + 1
        #     end_idx = min(i * chunk_size, n)
        #     #
        #     println(start_idx," ",end_idx)
        #     for x=start_idx:end_idx
        #         local M
        #         M = d2Uf[x]'*Va[x]*Uf[x] + 2*dUf[x]'*Va[x]*dUf[x] + Uf[x]'*Va[x]*d2Uf[x]
        #         for i=1:dim
        #             SmoothSum[i] += M[i,i]^2
        #         end
        #     end
        # end

        SmoothSum = zeros(Float64,dim)
        #
        for x=1:lastindex(r)
            #
            M = d2U[x]'*Va[x]*U[x] + 2*dU[x]'*Va[x]*dU[x] + U[x]'*Va[x]*d2U[x]
            #
            for i=1:dim
                for j=i:dim
                    SmoothSum[i] += M[i,j]^2
                end
            end
        end
        #
        return sqrt(sum(SmoothSum))
    end
    #
    ## switching generator matrices via sigmoids with r-dependent gammas
    function matMix_variable_gamma(x, A, B, r0_vec, g_morphed, idx)
        #
        ## initialise matrix C
        C = zeros(Float64,dim,dim)
        #
        ## assume row stacking indices
        k = 0
        #
        for i=1:dim
            for j=i+1:dim
                #
                k += 1
                #
                f = sigmoid_variable_gamma(x, g_morphed[k][idx], r0_vec[k])[1]
                #
                C[i,j] = f * B[i,j] + (1 - f) * A[i,j] 
                C[j,i] = f * B[j,i] + (1 - f) * A[j,i]
            end
        end
        #
        return C
    end
    #
    ## cost function for the morphing fit
    function cost_morph(rgrid, Va, Kf, Kb, dim, order, r0, gamma_refs, p)
        #
        gamma_par = [p[ (i-1) * (order+1) + 1: i * (order+1)] for i=1:Int(dim*(dim-1)/2)]
        #
        ## morph gammas
        k = 0
        g_morphed = []
        #
        g_morphed_boundary = []
        for i=1:dim
            for j=i+1:dim
                #
                k += 1
                #
                # gamma_par_ij = round.(gamma_par[k], sigdigits=2)
                #
                push!(g_morphed, morph_const_pd24(r, gamma_refs[k], r0[k], gamma_par[k]))
                #
                push!(g_morphed_boundary, morph_const_pd24(region, gamma_refs[k], r0[k], gamma_par[k]))
            end
        end
        #
        if any([any(g_morphed[i] .<= 0.0) for i=1:dim])
            return 1e100
        end
        #
        ## check if the sigmoids break boundary conditions
        switchMat_LBC = []
        switchMat_RBC = []
        k = 0
        for i=1:dim
            for j=i+1:dim
                #
                k += 1
                #
                fL = sigmoid(region[1], g_morphed_boundary[k][1], r0[k])
                fR = sigmoid(region[2], g_morphed_boundary[k][2], r0[k])
                #
                if any([fL > µ, fR < 1 - µ])
                    return 1e100
                end
            end
        end
        #
        ## compute mixing of Kf and Kb
        K  = map( x -> matMix_variable_gamma(rgrid[x], Kf[x], Kb[x], r0, g_morphed, x), collect(1:lastindex(rgrid)))
        #
        ## compute the corresponding unitary matrix
        U  = real(exp.(K))
        #
        ## compute smoothness of the diabatic curves induced by U (ignore adiabatic smoothness/derivatives)
        smoothness =  SmoothnessByU(rgrid,U,Va,dim)
        #
        return smoothness
    end
    #
    ## inversion of sigmoid for gamma
    function gamma_sigmoid(rref, f, r0)
        return -log((1/(f))-1)/(rref-r0)
    end
    #
    function sigmoid_square(r, rmin, rmax, rref; µ = 1e-2)
        #
        ##
        rLmid = rmin + 0.25*(rref - rmin)
        rRmid = rmax - 0.25*(rmax - rref)
        #
        ## find gamma to ensure the sigmoid is zero at boundaries
        gL = gamma_sigmoid(rmin, µ, rLmid)
        gR = gamma_sigmoid(rmax, 1-µ, rRmid)
        #
        g = maximum([gL, gR])
        #
        ## construct the sigmoid-square
        return sigmoid(r, gL, rLmid)*(1 - sigmoid(r, gR, rRmid))
    end
    #                  
    ################## INITIALISATION OF ROUTINE PARAMETERS ####################
    Nel = Int(dim*(dim-1)/2)
    #
    ## convert matrices to vector format
    Uf = [Uf[i,:,:] for i=1:lastindex(Uf[:,1,1])]
    Ub = [Ub[i,:,:] for i=1:lastindex(Ub[:,1,1])]
    Va = [Va[i,:,:] for i=1:lastindex(Va[:,1,1])]
    #
    ## compute generators for each solution
    Kf = real(log.(Uf))
    Kb = real(log.(Ub))
    #
    ## find average of peak positions (the elements are closest) to set the sigmoid centres
    ∂K = Kf .- Kb
    r0 = []
    r0_bounds = []
    g0 = []
    for i=1:dim
        for j=i+1:dim
            #
            ## compute position where f & b generator elements vary the quickest
            ∂Kij = [d[i,j] for d in ∂K]
            deriv_∂Kij  = FiniteDifference(rgrid, ∂Kij, 1)
            ij_peak_idx = argmax(abs.(deriv_∂Kij))
            push!(r0, rgrid[ij_peak_idx])
            #
            ## compute HWHM of this for fitting range
            deriv_∂Kij_peak   = maximum(abs.(deriv_∂Kij))
            # println(deriv_∂Kij_peak,deriv_∂Kij)
            plt.figure()
            plt.plot(r,deriv_∂Kij)
            #
            HWHM_end_idx   = ij_peak_idx + (findfirst(x -> x == 0, [abs(a) >= 0.5 * deriv_∂Kij_peak for a in deriv_∂Kij[ij_peak_idx:end]]) - 1) - 1
            HWHM_start_idx = length(reverse(∂Kij[1:ij_peak_idx])) - (findfirst(x -> x == 0, [abs(a) >= 0.5 * deriv_∂Kij_peak for a in reverse(deriv_∂Kij[1:ij_peak_idx]) ]) - 1) + 1
            #
            push!(r0_bounds, [rgrid[HWHM_start_idx], rgrid[HWHM_end_idx]])
            #
            ## now estimate mimimum value for gamma
            g_thresh_left  = [gamma_sigmoid(region[1],     µ, rgrid[HWHM_start_idx]), gamma_sigmoid(region[1],     µ, r0[end]), gamma_sigmoid(region[1],     µ, rgrid[HWHM_end_idx])]
            g_thresh_right = [gamma_sigmoid(region[2],   1-µ, rgrid[HWHM_start_idx]), gamma_sigmoid(region[2],   1-µ, r0[end]), gamma_sigmoid(region[2],   1-µ, rgrid[HWHM_end_idx])]
            #
            push!(g0, maximum([g_thresh_left..., g_thresh_right...]))
            #
        end
    end
    #
    ## initialise guess for the morphing parameters
    p_guess = []
    guess_ij = zeros(Float64, order + 1)
    guess_ij[1] = 1.0
    #
    for i=1:Nel
        #
        push!(p_guess, 1.0)
        #
        for j=1:order
            push!(p_guess, 0.0)
        end
    end
    #
    g_guess = [g0...,r0...]
    #
    ## minimiser options 
    options = Optim.Options(
    x_tol = 1e-4,
    show_trace = true
    )
    #                  
    ############ FITTING THE INITIAL SIGMOID SWITCHING FUNCTIONS ###############
    #
    # println(r0,"   ",Nel,"  ",[g0...,r0...])
    o_ = optimize(p -> cost_init_gamma(rgrid, Va, Kf, Kb, dim, p[Nel+1:end], p[1:Nel], g0, r0_bounds), [g0...,r0...], options)
    sig_refs = Optim.minimizer(o_)
    gamma_refs = sig_refs[1:Nel]
    r0 = sig_refs[Nel+1:end]
    # gamma_refs = [7.985958500435945, 5.348656061575216, 5.218273496037558]
    # r0 = [1.212984790628298,1.239056969450167,1.2157919237239567]
    #
    println()
    println("__OPTIMISED REFERENCE GAMMA PARAMETERS__")
    k = 0
    for i=1:dim
        for j=i+1:dim
            k += 1
            println("γ_",i,j," = ", gamma_refs[k],"\tr0 = ",r0[k])
        end
    end
    println()
    #                 
    ############## FITTING THE MORPHED SIGMOID GAMMA PARAMETERS ################
    #
    o_ = optimize(p -> cost_morph(rgrid, Va, Kf, Kb, dim, order, r0, gamma_refs, p), [p_guess...], options)
    gamma_pars = Optim.minimizer(o_)
    # println("COST = ", cost_morph(rgrid, Va, Kf, Kb, dim, order, r0, gamma_refs, gamma_pars))
    #
    gamma_pars = [gamma_pars[ (i-1) * (order+1) + 1: i * (order+1)] for i=1:Nel]
    # 
    # gamma_pars = [[9.761812719304206,-0.0691116594691795,-0.32582505708266024],[0.058901506522896004,-2.384240857909938,0.3663287614513354],[1.1768121977434012,1.1937730008615215,1.843439941798288]]
    println()
    println("__OPTIMISED GAMMA MORPHING PARAMETERS__")
    println()
    k = 0
    for i=1:dim
        for j=i+1:dim
            k += 1
            println("γ_",i,j,":")
            println("r0 = ", r0[k])
            println("γ0 = ", gamma_refs[k])
            for d=1:order+1
                println("B",d-1," = ", gamma_pars[k][d])
            end
            println()
        end
    end
    println()
    #
    rsolve = DENSE_GRID_OBJECTS[1]
    Uf_dense_grid = DENSE_GRID_OBJECTS[2]
    Ub_dense_grid = DENSE_GRID_OBJECTS[3]
    Uf_dense_grid = [ Uf_dense_grid[i,:,:] for i=1:lastindex(Uf_dense_grid[:,1,1])]
    Ub_dense_grid = [ Ub_dense_grid[i,:,:] for i=1:lastindex(Ub_dense_grid[:,1,1])]

    Kf_dense_grid = real(log.(Uf_dense_grid))
    Kb_dense_grid = real(log.(Ub_dense_grid))
    #
    ## compute morphed gammass
    k = 0
    g_morphed = []
    #
    g_morphed_dense_grid = []
    #
    for i=1:dim
        for j=i+1:dim
            #
            k += 1
            #
            # gamma_par_ij = round.(gamma_par[k], sigdigits=2)
            #
            println(gamma_refs[k], r0[k])
            push!(g_morphed, morph_const_pd24(rgrid, gamma_refs[k], r0[k], gamma_pars[k]))
            push!(g_morphed_dense_grid, morph_const_pd24(rsolve, gamma_refs[k], r0[k], gamma_pars[k]))
            #
        end
    end
    #
    ## compute the fitted generator matrix
    K  = map( x -> matMix_variable_gamma(rgrid[x], Kf[x], Kb[x], r0, g_morphed, x), collect(1:lastindex(rgrid)))
    #
    ## compute the corresponding unitary matrix
    U  = real(exp.(K))
    #
    Vd = adjoint.(U) .* Va .* U
    #
    dU = FiniteDiff_MatDerivative(rgrid,U,dim,1) # computes derivative of matrix
    UdU = map(i -> U[i]*dU[i]',collect(1:size(rgrid)[1]))
    #
    ## compute matrices on the dense grid
    K_dense_grid  = map( x -> matMix_variable_gamma(rsolve[x], Kf_dense_grid[x], Kb_dense_grid[x], r0, g_morphed_dense_grid, x), collect(1:lastindex(rsolve)))
    #
    U_dense_grid  = real(exp.(K_dense_grid))
    #
    dU_dense_grid = FiniteDiff_MatDerivative(rsolve,U_dense_grid,dim,1) # computes derivative of matrix
    #
    W_regularised = map(i -> U_dense_grid[i]*dU_dense_grid[i]',collect(1:size(rsolve)[1]))
    #
    return U, dU, UdU, K, Vd, W_regularised

end
#
function regularise(rgrid, Uf, Ub, Va, dim, order, region; DENSE_GRID_OBJECTS = false, µ = 0.1)
    rgrid = collect(rgrid)
    #                  
    ######################### FUNCTIONS FOR SUBROUTINE #########################
    #
    ## switching generator matrices with #dim sigmoids
    function matMix_const_gamma(x, A, B, r0_vec, gamma_vec)
        #
        ## initialise matrix C
        C = zeros(Float64,dim,dim)
        #
        ## assume row stacking indices
        k = 0
        #
        for i=1:dim
            for j=i+1:dim
                #
                k += 1
                #
                f = sigmoid(x, gamma_vec[k], r0_vec[k])
                #
                C[i,j] = f * B[i,j] + (1 - f) * A[i,j] 
                C[j,i] = f * B[j,i] + (1 - f) * A[j,i]
            end
        end
        #
        return C
    end
    #
    ## sigmoid with r-dependent gamma
    function sigmoid_variable_gamma(r,gamma,r0)
        return map(x -> 1/(1+exp(-gamma[x]*(r[x]-r0))), collect(1:lastindex(r)))
    end
    #
    ## morphing function
    function polynom_decay_24(r, re, B; Binf = 1, beta2 = 0.1, beta4 = 0.02, p = 4)
        #
        ## Surkus variable
        y = (r^p-re^p)/(r^p+re^p)
        #
        ## reduced variable
        z = (r-re) * exp(-beta2*(r-re)^2-beta4*(r-re)^4) # sigmoid_square(r, region[1], region[2], re)
        #
        expansion = sum([B[k] * z^(k-1) * (1-y) for k=1:size(B)[1]]) + Binf * y
        #
        return expansion
    end
    #
    ## morph constant to a curve via multiplication with f =  polynom_decay_24
    function morph_const_pd24(r, fconst, re, B; Binf = 1, beta2 = 0.1, beta4 = 0.02, p = 4)
        #
        return fconst .* map(x -> polynom_decay_24(x, re, B, Binf = Binf, beta2 = beta2, beta4 = beta4, p = p), r)
    end
    #
    ## finite difference (high order accuracy) for matrices
    function FiniteDiff_MatDerivative(x, M, dim, d_order)
        dM = zeros(Float64,length(x),dim,dim)
        #
        ## extract elements of U, spline them, form splined matrix and derivative
        for i=1:dim
            for j=1:dim
                #
                Mij = [u[i,j] for u in M]
                #
                dMij = FiniteDifference(x,Mij,d_order)
                #
                dM[:,i,j] = dMij
            end
        end
        #
        return [dM[idx,:,:] for idx=1:lastindex(x)]
    end
    #
    ## smoothness of the dynamics induced into diabats by U (omitting adiabats)
    function SmoothnessByU(r,U,Va,dim)
        #
        # if Threads.nthreads() >=2
        #     #
        #     derivs = Vector{Vector{SMatrix{dim, dim, Float64}}}(undef, 2)
        #     #
        #     @threads for i in 1:2
        #         derivs[i] = Vector{SMatrix{dim, dim, Float64}}(FiniteDiff_MatDerivative(r,U,dim,i))
        #     end
        #     #
        #     dU, d2U = derivs[1], derivs[2]
        # else
        dU  = Vector{SMatrix{dim, dim, Float64}}(FiniteDiff_MatDerivative(r,U,dim,1))
        d2U = Vector{SMatrix{dim, dim, Float64}}(FiniteDiff_MatDerivative(r,U,dim,2))
        # end
        #
        U  = Vector{SMatrix{dim, dim, Float64}}(U)
        Va = Vector{SMatrix{dim, dim, Float64}}(Va)
        #
        # n = length(r)
        # num_threads = min(Threads.nthreads(), n)
        # chunk_size = ceil(Int, n / num_threads)
        
        # SmoothSum = zeros(Float64,dim)
    
        # @threads for i in 1:num_threads
        #     start_idx = (i - 1) * chunk_size + 1
        #     end_idx = min(i * chunk_size, n)
        #     #
        #     println(start_idx," ",end_idx)
        #     for x=start_idx:end_idx
        #         local M
        #         M = d2Uf[x]'*Va[x]*Uf[x] + 2*dUf[x]'*Va[x]*dUf[x] + Uf[x]'*Va[x]*d2Uf[x]
        #         for i=1:dim
        #             SmoothSum[i] += M[i,i]^2
        #         end
        #     end
        # end

        SmoothSum = zeros(Float64,dim)
        #
        for x=1:lastindex(r)
            #
            M = d2U[x]'*Va[x]*U[x] + 2*dU[x]'*Va[x]*dU[x] + U[x]'*Va[x]*d2U[x]
            #
            for i=1:dim
                for j=i:dim
                    SmoothSum[i] += M[i,j]^2
                end
            end
        end
        #
        return sqrt(sum(SmoothSum))
    end
    #
    ## switching generator matrices via sigmoids with r-dependent gammas
    function matMix_variable_gamma(x, A, B, r0_vec, g_morphed, idx)
        #
        ## initialise matrix C
        C = zeros(Float64,dim,dim)
        #
        ## assume row stacking indices
        k = 0
        #
        for i=1:dim
            for j=i+1:dim
                #
                k += 1
                #
                f = sigmoid_variable_gamma(x, g_morphed[k][idx], r0_vec[k])[1]
                #
                C[i,j] = f * B[i,j] + (1 - f) * A[i,j] 
                C[j,i] = f * B[j,i] + (1 - f) * A[j,i]
            end
        end
        #
        return C
    end
    #
    ## slices array of switching function parameters into subsequent blocks based on number of parameters
    function parameter_idx_splicer(Nparams,idx)
        #
        if idx == 1
            return 1, 1 + (Nparams[1] - 1)
        else
            idx_start = sum(Nparams[1:idx-1])+1
            #
            idx_end   = idx_start + (Nparams[idx] - 1)
            #
            return idx_start, idx_end
        end
    end
    #
    ## cost function for fitting switching function
    function cost(rgrid, Va, Kf, Kb, dim, bounds, fit_flags, pidx, p)
        #
        ## extract parameters, check if to be fit, check for bounds, and compute morphing function
        row_col_counter = 0
        #
        r0 = []
        #
        g_morphed = []
        #
        g_morphed_boundary = []
        for i=1:dim
            for j=i+1:dim
                #
                ## extract parameter for element Fij
                row_col_counter += 1
                #
                idx_start, idx_end = pidx[row_col_counter]
                #
                param_ij = p[idx_start:idx_end]
                #
                ## check to see if any parameters fall outside of the boundary windows, if so return huge cost
                parameter_bounds = bounds[row_col_counter] #SwitchingFunction[[i,j]].bounds
                outOfBounds  = [(param_ij[i] < parameter_bounds[i][1])|(param_ij[i] > parameter_bounds[i][2]) for i=1:lastindex(param_ij)]
                #
                if any(outOfBounds .== 1)
                    return 1e100
                end
                #
                ## check which parameters are fit, if they are not to be fit set them equal to their user defined value
                fit_flag = fit_flags[row_col_counter] #SwitchingFunction[[i,j]].fit
                #
                for (idx, flag) in enumerate(fit_flag)
                    if flag != 1
                        param_ij[idx] = SwitchingFunction[[i,j]].Rval[idx]
                    end
                end
                #
                ## assign parameters
                g0_ij    = param_ij[1]
                r0_ij    = param_ij[2]
                p_ij     = param_ij[3]
                beta2_ij = param_ij[4]
                beta4_ij = param_ij[5]
                B_ij     = param_ij[6:end]
                #
                ## compute the morphed sigmoid width parameter
                push!(g_morphed, morph_const_pd24(r, g0_ij, r0_ij, B_ij, beta2 = beta2_ij, beta4 = beta4_ij, p = p_ij))
                #
                push!(g_morphed_boundary, morph_const_pd24(region, g0_ij, r0_ij, B_ij, beta2 = beta2_ij, beta4 = beta4_ij, p = p_ij))
                #
                push!(r0,r0_ij)
            end
        end
        #
        ## check if the sigmoids break boundary conditions
        row_col_counter = 0
        for i=1:dim
            for j=i+1:dim
                #
                row_col_counter += 1
                #
                fL = sigmoid(region[1], g_morphed_boundary[row_col_counter][1], r0[row_col_counter])
                fR = sigmoid(region[2], g_morphed_boundary[row_col_counter][2], r0[row_col_counter])
                #
                if any([fL > µ, fR < 1 - µ])
                    return 1e100
                end
            end
        end
        #
        ## compute mixing of Kf and Kb
        K  = map( x -> matMix_variable_gamma(rgrid[x], Kf[x], Kb[x], r0, g_morphed, x), collect(1:lastindex(rgrid)))
        #
        ## compute the corresponding unitary matrix
        U  = real(exp.(K))
        #
        ## compute smoothness of the diabatic curves induced by U (ignore adiabatic smoothness/derivatives)
        smoothness =  SmoothnessByU(rgrid,U,Va,dim)
        #
        return smoothness
    end
    #
    ## inversion of sigmoid for gamma
    function gamma_sigmoid(rref, f, r0)
        return -log((1/(f))-1)/(rref-r0)
    end
    #
    function sigmoid_square(r, rmin, rmax, rref; µ = 1e-2)
        #
        ##
        rLmid = rmin + 0.25*(rref - rmin)
        rRmid = rmax - 0.25*(rmax - rref)
        #
        ## find gamma to ensure the sigmoid is zero at boundaries
        gL = gamma_sigmoid(rmin, µ, rLmid)
        gR = gamma_sigmoid(rmax, 1-µ, rRmid)
        #
        g = maximum([gL, gR])
        #
        ## construct the sigmoid-square
        return sigmoid(r, gL, rLmid)*(1 - sigmoid(r, gR, rRmid))
    end
    #
    function initialise_switches_black_box(∂K, rgrid)
        #
        ## find average of peak positions (the elements are closest) to set the sigmoid centres
        r0 = []
        r0_bounds = []
        g0 = []
        for i=1:dim
            for j=i+1:dim
                if [i,j] ∉ keys(SwitchingFunction)
                    #
                    ## compute position where f & b generator elements vary the quickest
                    ∂Kij = [d[i,j] for d in ∂K]
                    deriv_∂Kij  = FiniteDifference(rgrid, ∂Kij, 1)
                    ij_peak_idx = argmax(abs.(deriv_∂Kij))
                    push!(r0, rgrid[ij_peak_idx])
                    #
                    ## compute HWHM of this for fitting range
                    deriv_∂Kij_peak   = maximum(abs.(deriv_∂Kij))
                    # println(deriv_∂Kij_peak,deriv_∂Kij)
                    # plt.figure()
                    # plt.plot(r,deriv_∂Kij)
                    #
                    HWHM_end_idx   = ij_peak_idx + (findfirst(x -> x == 0, [abs(a) >= 0.5 * deriv_∂Kij_peak for a in deriv_∂Kij[ij_peak_idx:end]]) - 1) - 1
                    HWHM_start_idx = length(reverse(∂Kij[1:ij_peak_idx])) - (findfirst(x -> x == 0, [abs(a) >= 0.5 * deriv_∂Kij_peak for a in reverse(deriv_∂Kij[1:ij_peak_idx]) ]) - 1) + 1
                    #
                    push!(r0_bounds, [rgrid[HWHM_start_idx], rgrid[HWHM_end_idx]])
                    #
                    ## now estimate mimimum value for gamma
                    g_thresh_left  = [gamma_sigmoid(region[1],     µ, rgrid[HWHM_start_idx]), gamma_sigmoid(region[1],     µ, r0[end]), gamma_sigmoid(region[1],     µ, rgrid[HWHM_end_idx])]
                    g_thresh_right = [gamma_sigmoid(region[2],   1-µ, rgrid[HWHM_start_idx]), gamma_sigmoid(region[2],   1-µ, r0[end]), gamma_sigmoid(region[2],   1-µ, rgrid[HWHM_end_idx])]
                    #
                    push!(g0, maximum([g_thresh_left..., g_thresh_right...]))
                    #
                else
                    push!(r0, mean(region))
                    push!(r0_bounds, [-1e100,1e100])
                    push!(g0, 1e-16)
                end
            end
        end
        #
        g_guess = [g0...,r0...]
        #
        return g0, r0, r0_bounds, g_guess
    end
    ################## INITIALISATION OF ROUTINE PARAMETERS ####################
    Nel = Int(dim*(dim-1)/2)
    #
    ## convert matrices to vector format
    Uf = [Uf[i,:,:] for i=1:lastindex(Uf[:,1,1])]
    Ub = [Ub[i,:,:] for i=1:lastindex(Ub[:,1,1])]
    Va = [Va[i,:,:] for i=1:lastindex(Va[:,1,1])]
    #
    ## compute generators for each solution & their difference matrix
    Kf = real(log.(Uf))
    Kb = real(log.(Ub))
    ∂K = Kf .- Kb
    #
    ## initialise switching function parameters
    g0, r0, r0_bounds, g_guess = initialise_switches_black_box(∂K,rgrid) # black box initialisation
    #
    ## row/column counter, e.g. 1,2,3,... -> [1,2], [1,3], ... , [2,3],[2,4], ...
    row_col_counter = 0
    #
    Nparams = []
    params_total = []
    param_idx = []
    #
    ## set user defined parameters
    for i=1:dim
        for j=i+1:dim
            row_col_counter += 1
            if [i,j] in keys(SwitchingFunction)
                param_ij = SwitchingFunction[[i,j]].Rval
            else
                switch_ij = SWITCH([i,j],
                                    "switch",
                                    ["g0", "r0","p","beta2","beta4","B0","B1","B2"],
                                    [g0[row_col_counter], r0[row_col_counter], 4, 0.1, 0.02, 1, 0, 0],
                                    [1,1,0,0,0,1,1,1],
                                    [[g0[row_col_counter]*0.99,1e100],r0_bounds[row_col_counter],[-1e100,1e100],[-1e100,1e100],[-1e100,1e100],[-1e100,1e100][-1e100,1e100],[-1e100,1e100]],
                                    [g0[row_col_counter], r0[row_col_counter], 4, 0.1, 0.02, 1, 0, 0]
                                   )
                #
                SwitchingFunction[[i,j]] = switch_ij
                #
                param_ij = [g0[row_col_counter], r0[row_col_counter], 4, 0.1, 0.02, 1, 0, 0]
            end
            #
            push!(Nparams,length(param_ij))
            #
            append!(params_total, param_ij)
            #
            push!(param_idx,parameter_idx_splicer(Nparams,row_col_counter))
        end
    end
    #
    ## minimiser options 
    options = Optim.Options(
    x_tol = 1e-4,
    show_trace = true
    )
    #
    ##################### FITTING THE SWITCHING FUNCTIONS ######################
    #
    o_ = optimize(p -> cost(rgrid, 
                            Va, 
                            Kf, 
                            Kb,
                            dim, 
                            [SwitchingFunction[[i,j]].bounds for i=1:dim-1 for j=i+1:dim], 
                            [SwitchingFunction[[i,j]].fit    for i=1:dim-1 for j=i+1:dim], 
                            param_idx, 
                            p), [params_total...], options)
    popt = Optim.minimizer(o_)
    #
    ## set the fitted parameters in the switching function fields and print the fitted parameters
    row_col_counter = 0
    #
    println()
    println("__OPTIMISED REFERENCE GAMMA PARAMETERS__")
    for i=1:dim
        for j=i+1:dim
            row_col_counter += 1
            #
            idx_start, idx_end = param_idx[row_col_counter]
            #
            fitted_params_ij = popt[idx_start:idx_end]
            #
            SwitchingFunction[[i,j]].fitted_parameters = fitted_params_ij
            #
            ## print the parameters for the user
            println("< ",i," | F | ",j," > :")
            println("---------------")
            println("g0"   , i, j, " = ", fitted_params_ij[1])
            println("r0"   , i, j, " = ", fitted_params_ij[2])
            println("p"    , i, j, " = ", fitted_params_ij[3])
            println("beta2", i, j, " = ", fitted_params_ij[4])
            println("beta4", i, j, " = ", fitted_params_ij[5])
            for d=1:lastindex(fitted_params_ij[6:end])
                println("B", d-1, " = ", fitted_params_ij[6:end][d])
            end
            println()
        end
    end
    #
    ################## COMPUTING THE REGULARISING CORRECTION ###################
    #
    ## extract the AtDTs on the solution grid
    rsolve = DENSE_GRID_OBJECTS[1]
    Uf_dense_grid = DENSE_GRID_OBJECTS[2]
    Ub_dense_grid = DENSE_GRID_OBJECTS[3]
    Uf_dense_grid = [ Uf_dense_grid[i,:,:] for i=1:lastindex(Uf_dense_grid[:,1,1])]
    Ub_dense_grid = [ Ub_dense_grid[i,:,:] for i=1:lastindex(Ub_dense_grid[:,1,1])]
    #
    ## compute the generator matrices on the solution grid
    Kf_dense_grid = real(log.(Uf_dense_grid))
    Kb_dense_grid = real(log.(Ub_dense_grid))
    #
    ## compute morphed gammas
    g_morphed = []
    #
    g_morphed_dense_grid = []
    #
    r0 = []
    #
    for i=1:dim
        for j=i+1:dim
            #
            param_ij = SwitchingFunction[[i,j]].fitted_parameters
            #
            g0_ij    = param_ij[1]
            r0_ij    = param_ij[2]
            p_ij     = param_ij[3]
            beta2_ij = param_ij[4]
            beta4_ij = param_ij[5]
            B_ij     = param_ij[6:end]
            #
            ## compute the morphed sigmoid width parameter
            push!(g_morphed, morph_const_pd24(rgrid, g0_ij, r0_ij, B_ij, beta2 = beta2_ij, beta4 = beta4_ij, p = p_ij))
            #
            push!(g_morphed_dense_grid, morph_const_pd24(rsolve, g0_ij, r0_ij, B_ij, beta2 = beta2_ij, beta4 = beta4_ij, p = p_ij))
            #
            push!(r0, r0_ij)
        end
    end
    #
    ## compute the fitted generator matrix
    K  = map( x -> matMix_variable_gamma(rgrid[x], Kf[x], Kb[x], r0, g_morphed, x), collect(1:lastindex(rgrid)))
    #
    ## compute the corresponding unitary matrix
    U  = real(exp.(K))
    #
    Vd = adjoint.(U) .* Va .* U
    #
    dU = FiniteDiff_MatDerivative(rgrid,U,dim,1) # computes derivative of matrix
    UdU = map(i -> U[i]*dU[i]',collect(1:size(rgrid)[1]))
    #
    ## compute matrices on the dense grid
    K_dense_grid  = map( x -> matMix_variable_gamma(rsolve[x], Kf_dense_grid[x], Kb_dense_grid[x], r0, g_morphed_dense_grid, x), collect(1:lastindex(rsolve)))
    #
    U_dense_grid  = real(exp.(K_dense_grid))
    #
    dU_dense_grid = FiniteDiff_MatDerivative(rsolve,U_dense_grid,dim,1) # computes derivative of matrix
    #
    W_regularised = map(i -> U_dense_grid[i]*dU_dense_grid[i]',collect(1:size(rsolve)[1]))
    #
    return U, dU, UdU, K, Vd, W_regularised
end
#
function N_state_diabatisation(block)
    dim = length(block)
    block =  sort(block, by = t -> t[1])   
    renumbered_states =  Dict{Int, Int}()
    for i=1:dim
        renumbered_states[block[i]] = i
    end
    #
    ## compute number of grid splits based on user input on grid resolution
    nsplit = 0
    ∂r = 100
    spacing = []
    while ∂r > Calculation["method"].grid_resolution
        # global nsplit, ∂r
        r_ = LinRange(r[1], r[end], gridSplitter(size(r)[1],nsplit))
        ∂r = r_[2]-r_[1]
        push!(spacing,∂r)
        nsplit+=1
    end
    nsplit = argmin(spacing .- Calculation["method"].grid_resolution) -1
    #
    ## compute the extended grid
    rsolve, ∆ = smoothgrid(-1000,1000,r,nsplit,1000,1000)
    #
    ## initialse the evolving NAC matrix
    evoNACMat = zeros(length(rsolve),dim,dim)
    #
    for key in keys(NonAdiabaticCoupling)
        # if all(x -> x in block, key)
        if check_symmetry_and_states(Int(key[1]),Int(key[2]))
            local _, NAC
            #
            i, j = renumbered_states[key[1]], renumbered_states[key[2]]
            #
            _, NAC = ComputeProperty(NonAdiabaticCoupling[key], custom_grid = rsolve, evolution_grid = true)
            evoNACMat[:,i,j] =   NAC .* NonAdiabaticCoupling[key].factor
            evoNACMat[:,j,i] = - NAC .* NonAdiabaticCoupling[key].factor
        end
    end
    #
    ## now perform evolution
    if  occursin("forward",lowercase(Calculation["method"].diabatisation))
        U, dU, UdU = Forward_Evolution(rsolve, evoNACMat, dim, Calculation["method"].l_boundary_condition)
        #
        U  = SplineMat(rsolve,  U, r)
        dU = SplineMat(rsolve, dU, r)
    #
    elseif occursin("backward",lowercase(Calculation["method"].diabatisation))
        #
        ## determine the right boundary condition
        if Calculation["method"].r_boundary_condition == false
            Uf, _, _ = Forward_Evolution(rsolve, evoNACMat, dim, Calculation["method"].l_boundary_condition)
            #
            r_boundary_condition = find_closest_permutation(Uf[end,:,:], dim)
        else
            r_boundary_condition = Calculation["method"].r_boundary_condition
        end
        #
        U, dU, UdU = Backward_Evolution(rsolve, evoNACMat, dim, r_boundary_condition)
        #
        U  = SplineMat(rsolve,  U, r)
        dU = SplineMat(rsolve, dU, r)
    #
    elseif (Calculation["method"].regularisation != false)&(lowercase(Calculation["method"].diabatisation) == "evolution")
            #
            ## perform forward evolution
            Uf, dUf, UdUf = Forward_Evolution(rsolve, evoNACMat, dim, Calculation["method"].l_boundary_condition)
            #
            ## determine the right boundary condition
            if Calculation["method"].r_boundary_condition == false
                r_boundary_condition = find_closest_permutation(Uf[end,:,:], dim)
            else
                r_boundary_condition = Calculation["method"].r_boundary_condition
            end
            #
            ## perform backward evolution
            Ub, dUb, UdUb = Backward_Evolution(rsolve, evoNACMat, dim, r_boundary_condition)
            #
            ## determine the adiabatic property to regualrise to
            if Calculation["method"].regularisation == false
                adi_prop = Objects["potential"]
            else
                adi_prop = Objects[Calculation["method"].regularisation]
            end
            #
            ## initialse the adiabatic property matrix
            Pa = zeros(Float64,length(r),dim,dim)
            #
            for i in Calculation["method"].states
                for j in Calculation["method"].states
                    if check_symmetry_and_states(Int(i),Int(j))
                        row, col = renumbered_states[i], renumbered_states[j]
                        #
                        Pa[:,row,col] = adi_prop[:,i,j]
                    end
                end
            end
            #
            ## perform regularisation procedure
            U_, dU_, UdU_, K, Vd, W_regularised = regularise(r, 
                                                          SplineMat(rsolve, Uf, r), 
                                                          SplineMat(rsolve, Ub, r), 
                                                          Pa, 
                                                          dim, 
                                                          2, 
                                                          [r[1], r[end]], 
                                                          DENSE_GRID_OBJECTS = [rsolve, Uf, Ub])
            #
            UdU = W_regularised
            #
            U = zeros(length(r),dim,dim)
            dU = zeros(length(r),dim,dim)
            #
            for i in Calculation["method"].states
                for j in Calculation["method"].states
                    U[:,i,j] = [u[renumbered_states[i],renumbered_states[j]] for u in U_]
                    #
                    dU[:,i,j] = [x[renumbered_states[i],renumbered_states[j]] for x in dU_]
                end
            end
    end
    #
    return U, dU, UdU, rsolve, renumbered_states

end
#
function run_diabatiser(diabMethod)
    local U, dU, UdU, K_Matrix, diabatic_basis, Diabatic_Objects, psi_d, input_properties
    #
    # Calculation = diabatom["Calculation"]
    # Objects     = diabatom["Objects"]
    # Potential   = diabatom["Potential"]
    # Hamiltonian = diabatom["Hamiltonian"]
    #
    # diabMethod = lowercase(Calculation["method"].diabatisation)
    #
    dim = length(keys(Potential))
    # perform the diabatisation
    if (lowercase(diabMethod) == "2-state-approx")|(length(Calculation["method"].states)==2)
        U = fit_multi_diabat_2stateApprox(r, Objects["potential"])
        #
        UdU = Objects["nac"]
        #
        dU = map(i -> - UdU[i,:,:] * U[i,:,:]', collect(1:lastindex(r)))
    elseif  occursin("evolution",diabMethod)
        U_, dU, UdU_, rsolve, renumbered_states = N_state_diabatisation(Calculation["method"].states)
        #
        ## U_, dU_, has elements of U_[:,i,j], UdU_ has elements in list [x[i,j] for x in UdU_]
        #
        ## populate the full dimensional U and UdU matrix
        U = zeros(length(r),dim,dim)
        for i=1:dim
            U[:,i,i] .= 1.0
        end
        #
        UdU = zeros(length(rsolve),dim,dim)
        #
        for i in Calculation["method"].states
            for j in Calculation["method"].states
                U[:,i,j]   = U_[:,renumbered_states[i],renumbered_states[j]] # for u in U_]
                #
                if i!=j
                    UdU[:,i,j] = [w[renumbered_states[i],renumbered_states[j]] for w in UdU_]
                end
            end
        end
        #
        ## spline UdU onto the PEC grid
        UdU = SplineMat(rsolve, UdU, r)
        # #
        # ## create new regularised NAC object
        # Objects["regularised_nac"] = UdU
        # #
        # ## compute K matrix (<di/dr|dj/dr>)
        # K_Matrix_ = map(i -> -UdU[i,:,:] * UdU[i,:,:], collect(1:lastindex(r)))
        # K_Matrix = zeros(length(r),dim,dim)
        # for i=1:dim
        #     for j=1:dim
        #         K_Matrix[:,i,j] = [k[i,j] for k in K_Matrix_]
        #     end
        # end
        # Objects["K_matrix"] = K_Matrix
    end
    #
    ## create new regularised NAC object
    Objects["regularised_nac"] = UdU
    #
    ## compute K matrix (<di/dr|dj/dr>)
    K_Matrix_ = map(i -> -UdU[i,:,:] * UdU[i,:,:], collect(1:lastindex(r)))
    K_Matrix = zeros(length(r),dim,dim)
    for i=1:dim
        for j=1:dim
            K_Matrix[:,i,j] = [k[i,j] for k in K_Matrix_]
        end
    end
    Objects["K_matrix"] = K_Matrix
    #
    ## diabatic basis
    psi_d = zeros(Float64,dim)
    diabatic_basis = []
    for i=1:dim
        psi_d = zeros(Float64,dim)
        psi_d[i] = 1.0
        push!(diabatic_basis, map(x -> U[x,:,:]*psi_d, collect(1:length(r))))
    end
    #
    ## diabatic properties
    input_properties = unique([k[1] for k in keys(Hamiltonian)])
    #
    Diabatic_Objects = Dict()
    for p in input_properties
    if p != "NAC"
        #
        if p == "poten"
            #
            Adiabatic_PotMat = Objects["potential"]
            #
            Diabatic_PotMat = map(i -> U[i,:,:]' * Adiabatic_PotMat[i,:,:] * U[i,:,:], collect(1:lastindex(r)))
            #
            ##
            Diabatic_PropMat = zeros(Float64, length(r), dim, dim)
            #
            for i=1:dim
                for j=1:dim
                    Diabatic_PropMat[:,i,j] = [Vd[i,j] for Vd in Diabatic_PotMat]
                end
            end
            #
            Diabatic_Objects["potential"] = Diabatic_PropMat
        else
            Adiabatic_PropMat = Objects[lowercase(p)]
            #
            Diabatic_PropMat = zeros(Float64, length(r), dim, dim)
            #
            for i=1:dim
                for j=1:dim
                    Pd_ij =  Diabatic_Property_from_wavefunctions(diabatic_basis[i], diabatic_basis[j], Adiabatic_PropMat)
                    #
                    Diabatic_PropMat[:,i,j] = Pd_ij
                end
            end
            #
            Diabatic_Objects[p] = Diabatic_PropMat
        end
    end
    end
    #
    property_ordering = ["poten","nac","spin-orbit","lx","dipole"]
    #
    sorted_properties = [p for p in property_ordering if lowercase(p) in lowercase.(input_properties)]
    #
    return U, dU, UdU, K_Matrix, diabatic_basis, Diabatic_Objects, sorted_properties
end
#
function save_diabatisation(Objects, Diabatic_Objects, diabMethod, input_properties, fname; special_name = "", folder_path = false, intensity_flag = false, thresh_intes  = 1e-40, thresh_line   = 1e-40, thresh_bound  = 1e-4, thresh_dipole = 1e-8)
    function duo_calc_setup(io, atom1, atom2, dim, r, jrot, vmax)
        #
        ## rovibronic symmetry
        if atom1==atom2
            calc_symmetry = "C2v"
        elseif atom1!=atom2
            calc_symmetry = "Cs(M)"
        end
        #
        ## duo input setup
        println(io, "masses "*string(atom1)*" "*string(atom2))
        println(io, "")
        println(io, "nstates "*string(dim))
        println(io, "")
        println(io, "jrot "*string(jrot))
        println(io, "")
        println(io, "grid")
        println(io, "  npoints "*string(length(r)))
        println(io, "  range   "*string(r[1])*" "*string(r[end]))
        println(io, "  type 0")
        println(io, "end")
        println(io, "")
        println(io, "(PRINT_PECS_AND_COUPLINGS_TO_FILE)")
        println(io, "")
        println(io, "CONTRACTION")
        println(io, " vib")
        println(io, " vmax "*string(vmax))
        println(io, "END")
        println(io, "")
        println(io, "DIAGONALIZER")
        println(io, "SYEVR")
        println(io, "nroots "*Calculation["save"].nroots)
        println(io, "end")
        println(io, "")
        println(io, "symmetry "*string(calc_symmetry))
        println(io, "")
    end
    #
    function duo_intensity_block(io, fname, rep, jrot, Emax; intensity_flag = false, thresh_intes  = 1e-40, thresh_line   = 1e-40, thresh_bound  = 1e-4, thresh_dipole = 1e-8)
        if intensity_flag == false
            println(io,"INTENSITY off")
        elseif intensity_flag == true
            println(io,"INTENSITY")
        end
        println(io,"absorption")
        println(io,"bound ")
        println(io,"THRESH_INTES  "*string(thresh_intes ))
        println(io,"THRESH_LINE   "*string(thresh_line  ))
        println(io,"thresh_bound  "*string(thresh_bound ))
        println(io,"thresh_dipole "*string(thresh_dipole))
        println(io,"TEMPERATURE   3000.0  (  RENORMALIZE)")
        println(io,"linelist  "*fname[1:end-4]*"_"*special_name*"_"*rep)
        println(io,"J, "*string(jrot))
        println(io,"freq-window    0, "*string(Emax))
        println(io,"energy low   -0.001, "*string(Emax)*", upper   -0.00, "*string(Emax))
        println(io,"")
    end
    #
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
    #
    function write_PEC(io,r, state, Potential_Matrix)
        V = Hamiltonian[("poten",state)]
        #
        name = V.name
        #
        symmetry = V.symmetry
        #
        L = V.lambda
        #
        S = V.mult
        #
        println(io, "poten " * string(state))
        println(io, "name '" * name * "'")
        println(io, "symmetry " * symmetry)
        println(io, "lambda " * string(L))
        println(io, "mult " * string(S))
        println(io, "type  grid")
        println(io, "values")
        #
        for idx=1:size(r)[1]
            x = r[idx]
            y = Potential_Matrix[idx, state, state]
            println(io, "\t$x\t$y")
        end
        #
        println(io, "end")
        println(io, "")
    end
    #
    function write_NAC(io,r, i, j, NAC_Matrix)
        W = Hamiltonian[("NAC",[i,j])]
        #
        name = W.name
        #
        V = Hamiltonian[("poten",i)]
        #
        symmetry = V.symmetry
        #
        L = V.lambda
        #
        S = V.mult
        #
        println(io, "NAC " * string(i) * " " * string(j))
        println(io, "name '" * name * "'")
        println(io, "symmetry " * symmetry)
        println(io, "lambda " * string(L))
        println(io, "mult " * string(S))
        println(io, "type  grid")
        println(io, "values")
        #
        for idx=1:size(r)[1]
            x = r[idx]
            y = NAC_Matrix[idx, i, j]
            println(io, "\t$x\t$y")
        end
        #
        println(io, "end")
        println(io, "")
    end
    #
    function write_W2(io,r,i,j,K_Matrix,factor)
        V = Hamiltonian[("poten",i)]
        #
        symmetry = V.symmetry
        #
        L = V.lambda
        #
        S = V.mult
        #
        name = "< "*string(i)*" | ∇ | "*string(j)*" >"
        #
        println(io, "diabatic " * string(i) * " " * string(j))
        println(io, "name '" * name * "'")
        println(io, "symmetry " * symmetry)
        println(io, "lambda " * string(L))
        println(io, "mult " * string(S))
        println(io, "type  grid")
        println(io, "factor "*string(factor))
        println(io, "values")
        #
        for idx=1:size(r)[1]
            x = r[idx]
            y = K_Matrix[idx, i, j]
            println(io, "\t$x\t$y")
        end
        #
        println(io, "end")
        println(io, "")
    end
    #
    function write_DC(io,r,i,j,Vd_Matrix)
        V = Hamiltonian[("poten",i)]
        #
        symmetry = V.symmetry
        #
        L = V.lambda
        #
        S = V.mult
        #
        name = "< "*string(i)*" | DC | "*string(j)*" >"
        #
        println(io, "diabat " * string(i) * " " * string(j))
        println(io, "name '" * name * "'")
        println(io, "symmetry " * symmetry)
        println(io, "lambda " * string(L))
        println(io, "mult " * string(S))
        println(io, "type  grid")
        println(io, "values")
        #
        for idx=1:size(r)[1]
            x = r[idx]
            y = Vd_Matrix[idx, i, j]
            println(io, "\t$x\t$y")
        end
        #
        println(io, "end")
        println(io, "")
    end
    #
    function write_SOC(io,r,i,j,SO_Matrix)
        SO = Hamiltonian[("spin-orbit",[i,j])]
        #
        Lz = SO.Lz
        Lz = [string(Lz[1])," ",string(Lz[2])]
        #
        Vi = Hamiltonian[("poten",i)]
        #
        symmetryi = Vi.symmetry
        #
        Li = Vi.lambda
        #
        Si = Vi.mult
        #
        Vj = Hamiltonian[("poten",j)]
        #
        symmetryj = Vj.symmetry
        #
        Lj = Vj.lambda
        #
        Sj = Vj.mult
        #
        symmetry = [symmetryi," ",symmetryj]
        #
        L = [string(Li)," ",string(Lj)]
        #
        S = [string(Si)," ",string(Sj)]
        #
        ## object name
        name = "< "*string(i)*" | SO | "*string(j)*" >"
        #
        if i != j
            println(io, "spin-orbit-x " * string(i) * " " * string(j))
        else
            println(io, "spin-orbit " * string(i) * " " * string(j))
        end
        #
        println(io, "name '" * name * "'")
        println(io, "symmetry " * join(symmetry))
        println(io, "lambda " * join(L))
        println(io, "mult " * join(S))
        println(io, "<x|Lz|y> " * join(Lz))
        println(io, "type  grid")
        println(io, "factor 1")
        println(io, "values")
        #
        for idx=1:size(r)[1]
            x = r[idx]
            y = SO_Matrix[idx, i, j]
            println(io, "\t$x\t$y")
        end
        #
        println(io, "end")
        println(io, "")
    end
    #
    function write_EAM(io,r,i,j,Lx_Matrix)
        Lx = Hamiltonian[("lx",[i,j])]
        #
        Lz = Lx.Lz
        Lz = [string(Lz[1])," ",string(Lz[2])]
        #
        Vi = Hamiltonian[("poten",i)]
        #
        # symmetryi = Vi.symmetry
        #
        Li = Vi.lambda
        #
        Si = Vi.mult
        #
        Vj = Hamiltonian[("poten",j)]
        #
        # symmetryj = Vj.symmetry
        #
        Lj = Vj.lambda
        #
        Sj = Vj.mult
        #
        # symmetry = [symmetryi," ",symmetryj]
        #
        L = [string(Li)," ",string(Lj)]
        #
        S = [string(Si)," ",string(Sj)]
        #
        ## object name
        name = "< "*string(i)*" | Lx | "*string(j)*" >"
        #
        println(io, "Lx " * string(i) * " " * string(j))
        println(io, "name '" * name * "'")
        # println(io, "symmetry " * join(symmetry))
        println(io, "lambda " * join(L))
        println(io, "mult " * join(S))
        println(io, "<x|Lz|y> " * join(Lz))
        println(io, "type  grid")
        println(io, "factor 1")
        println(io, "values")
        #
        for idx=1:size(r)[1]
            x = r[idx]
            y = Lx_Matrix[idx, i, j]
            println(io, "\t$x\t$y")
        end
        #
        println(io, "end")
        println(io, "")
    end
    #
    function write_dipole(io,r,i,j,DM_Matrix)
        dm = Hamiltonian[("dipole",[i,j])]
        #
        Lz = dm.Lz
        #
        Vi = Hamiltonian[("poten",i)]
        #
        # symmetryi = Vi.symmetry
        #
        Li = Vi.lambda
        #
        Si = Vi.mult
        #
        Vj = Hamiltonian[("poten",j)]
        #
        # symmetryj = Vj.symmetry
        #
        Lj = Vj.lambda
        #
        Sj = Vj.mult
        #
        # symmetry = [symmetryi," ",symmetryj]
        #
        L = [string(Li)," ",string(Lj)]
        #
        S = [string(Si)," ",string(Sj)]
        #
        Lz_arr = [string(Int(Li))*"i",string(Int(Lj))*"i"]
        if Li == 0.0
            Lz_arr[1] = "0"
        elseif Lj == 0
            Lz_arr[3] = "0"
        end
        #
        if Lz[1] == "N/A"
            Lz_arr[1] = string(Li)*"i"
        else
            Lz_arr[1] = Lz[1]
        end
        if Lz[2] == "N/A"
            Lz_arr[2] = string(Lj)*"i"
        else
            Lz_arr[2] = Lz[2]
        end
        #
        Lz = [string(Lz_arr[1])," ",string(Lz_arr[2])]
        #
        ## object name
        name = "< "*string(i)*" | DM | "*string(j)*" >"
        #
        if i != j
            println(io, "dipole-x " * string(i) * " " * string(j))
        else
            println(io, "dipole " * string(i) * " " * string(j))
        end
        #
        println(io, "name '" * name * "'")
        # println(io, "symmetry " * join(symmetry))
        println(io, "lambda " * join(L))
        println(io, "mult " * join(S))
        println(io, "<x|Lz|y> " * join(Lz))
        println(io, "type  grid")
        println(io, "factor 1")
        println(io, "values")
        #
        for idx=1:size(r)[1]
            x = r[idx]
            y = DM_Matrix[idx, i, j]
            println(io, "\t$x\t$y")
        end
        #
        println(io, "end")
        println(io, "")
    end
    #
    ## check if folder exists
    # if folder_path != false
    #     if isdir(folder_path)
    #         #
    #         ## current date and time
    #         current_datetime = now()
    #         #
    #         folder_path = folder_path*"_"*string(current_datetime)
    #         #
    #         mkdir(folder_path)
    #     else
    #         mkdir(folder_path)
    #     end
    # else
    #     folder_path = ""
    # end
    #
    ## specify atomic masses
    atom1, atom2 = Calculation["method"].atoms
    kinetic_factor =  KE_factor(atom1, atom2)
    #
    if lowercase(Calculation["save"].as) == "duo"
        #
        ## initialise quantities
        atom1, atom2 = Calculation["method"].atoms
        jrot = Calculation["save"].jrot
        vmax = Calculation["save"].vmax
        #
        kinetic_factor =  KE_factor(atom1, atom2)
        #
        ##
        for rep in ["adi","dia"]
            open(fname[1:end-4]*"_duo_"*special_name*"_"*rep*".inp", "w") do io
            #
            ## write the duo calculation setup
            duo_calc_setup(io, atom1, atom2, dim, r, jrot, vmax)
            #
            ## now write curves
            for property in input_properties
                if lowercase(property) == "poten"
                    for i=1:dim
                        if rep == "adi"
                            write_PEC(io,r, i, Objects["potential"])
                        else
                            write_PEC(io,r, i, Diabatic_Objects["potential"])
                            #
                            for j=i+1:dim
                                if (i in Calculation["method"].states)&(j in Calculation["method"].states)
                                    write_DC(io,r, i, j, Diabatic_Objects["potential"])
                                end
                            end
                        end
                    end
                    #
                    ## write diabatic couplings
                    if rep == "dia"
                        for i=1:dim
                            for j=i+1:dim
                                if (i in Calculation["method"].states)&(j in Calculation["method"].states)
                                    write_DC(io,r, i, j, Diabatic_Objects["potential"])
                                end
                            end
                        end
                    end
                #
                elseif (lowercase(property) == "nac")&(rep == "adi")
                    for i=1:dim
                        for j=i:dim
                                if i != j
                                    if ("NAC",[i,j]) in keys(Hamiltonian) #(i in Calculation["method"].states)&(j in Calculation["method"].states)
                                        if diabMethod == "evolution"
                                            write_NAC(io,r, i, j, Objects["regularised_nac"])
                                        else
                                            write_NAC(io,r, i, j, NACMat)
                                        end
                                    end
                                end
                                #
                                if (i in Calculation["method"].states)&(j in Calculation["method"].states)
                                    write_W2(io,r, i, j, Objects["K_matrix"], kinetic_factor)
                                end
                        end
                    end
                #
                elseif lowercase(property) == "spin-orbit"
                    for i=1:dim
                        for j=i:dim
                            if ("spin-orbit",[i,j]) in keys(Hamiltonian)
                                if rep == "adi"
                                    write_SOC(io,r, i, j, Objects["spin-orbit"])
                                else
                                    write_SOC(io,r, i, j, Diabatic_Objects["spin-orbit"])
                                end
                            end
                        end
                    end
                #
                elseif lowercase(property) == "dipole"
                    for i=1:dim
                        for j=i:dim
                            if ("dipole",[i,j]) in keys(Hamiltonian)
                                if rep == "adi"
                                    write_dipole(io,r, i, j, Objects["dipole"])
                                else
                                    write_dipole(io,r, i, j, Diabatic_Objects["dipole"])
                                end
                            end
                        end
                    end
                #
                elseif lowercase(property) == "lx"
                    for i=1:dim
                        for j=i:dim
                            if ("Lx",[i,j]) in keys(Hamiltonian)
                                if rep == "adi"
                                    write_EAM(io,r, i, j, Objects["lx"])
                                else
                                    write_EAM(io,r, i, j, Diabatic_Objects["lx"])
                                end
                            end
                        end
                    end
                end
            end
            #
            ## write intensity block
            Emax = maximum(Objects["potential"][end,:,:]) - minimum(Objects["potential"][:,1,1])
            duo_intensity_block(io,fname, rep, jrot, Emax, intensity_flag = intensity_flag, thresh_intes  = thresh_intes, thresh_line   = thresh_line, thresh_bound  = thresh_bound, thresh_dipole = thresh_dipole)
            #
            ##
            println(io,"END")
            println(io,"")
            println(io,"")
            println(io,"")
            println(io,"")
            println(io,"")
            end
        end
    end
    #
    if lowercase(Calculation["save"].as) == "dataframe"
        for rep in ["adi","dia"]
            #
            ## create DataFrames
            for property in input_properties
                if property == "poten"
                    df_pot = DataFrame()
                    df_pot[!,"R"] = r
                    #
                    df_DC = DataFrame()
                    df_DC[!,"R"] = r
                    #
                    Adiabatic_PotMat = Objects["potential"]
                    #
                    for i=1:dim
                        if rep == "adi"
                            df_pot[!,"V"*string(i)] = Adiabatic_PotMat[:,i,i]
                        else
                            df_pot[!,"V"*string(i)] = Diabatic_Objects["potential"][:,i,i]
                            #
                            for j=i+1:dim
                                if (i in Calculation["method"].states)&(j in Calculation["method"].states)
                                    df_DC[!,"<"*string(i)*"|Vd|"*string(j)*">"] = Diabatic_Objects["potential"][:,i,j]
                                end
                            end
                        end
                    end
                    #
                    CSV.write(rep*"_poten.dat", df_pot)
                    CSV.write("DC.dat", df_DC)
                #
                elseif (property == "nac")&(rep == "adi")
                    df_W = DataFrame()
                    df_W[!,"R"] = r
                    for i=1:dim
                        for j=i:dim
                                if i != j
                                    if ("NAC",[i,j]) in keys(Hamiltonian) #(i in Calculation["method"].states)&(j in Calculation["method"].states)
                                        if diabMethod == "evolution"
                                            df_W[!,"<"*string(i)*"|d/dr|"*string(j)*">"] = Objects["regularised_nac"][:,i,j]
                                        else
                                            df_W[!,"<"*string(i)*"|d/dr|"*string(j)*">"] = NACMat[:,i,j]
                                        end
                                    end
                                end
                                #
                                if (i in Calculation["method"].states)&(j in Calculation["method"].states)
                                    df_W[!,"<dΨ"*string(i)*"/dr|dΨ"*string(j)*"/dr>"] = Objects["K_matrix"][:,i,j] .* kinetic_factor
                                end
                        end
                    end
                    #
                    CSV.write("DDR.dat", df_W)
                #
                elseif property == "spin-orbit"
                    df_SO = DataFrame()
                    df_SO[!,"R"] = r
                    for i=1:dim
                        for j=i:dim
                            if ("spin-orbit",[i,j]) in keys(Hamiltonian)
                                if rep == "adi"
                                    df_W[!,"<"*string(i)*"|SO|"*string(j)*">"] = Objects["spin-orbit"][:,i,j]
                                else
                                    df_W[!,"<"*string(i)*"|SO|"*string(j)*">"] = Diabatic_Objects["spin-orbit"][:,i,j]
                                end
                            end
                        end
                    end
                    #
                    CSV.write("spin-orbit"*rep*".dat", df_SO)
                #
                elseif property == "dipole"
                    df_DM = DataFrame()
                    df_DM[!,"R"] = r
                    for i=1:dim
                        for j=i:dim
                            if ("dipole",[i,j]) in keys(Hamiltonian)
                                if rep == "adi"
                                    df_DM[!,"<"*string(i)*"|DM|"*string(j)*">"] = Objects["dipole"][:,i,j]
                                else
                                    df_DM[!,"<"*string(i)*"|DM|"*string(j)*">"] = Diabatic_Objects["dipole"][:,i,j]
                                end
                            end
                        end
                    end
                    #
                    CSV.write("dipole_"*rep*".dat", df_DM)
                #
                elseif property == "lx"
                    df_LX = DataFrame()
                    df_LX[!,"R"] = r
                    for i=1:dim
                        for j=i:dim
                            if ("Lx",[i,j]) in keys(Hamiltonian)
                                if rep == "adi"
                                    df_XLX[!,"<"*string(i)*"|Lx|"*string(j)*">"] = Objects["lx"][:,i,j]
                                else
                                    df_LX[!,"<"*string(i)*"|Lx|"*string(j)*">"] = Diabatic_Objects["lx"][:,i,j]
                                end
                            end
                        end
                    end
                    #
                    CSV.write("Lx"*rep*".dat", df_LX)
                end
            end
        end
    end
end