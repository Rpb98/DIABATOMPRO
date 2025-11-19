#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CALCULUS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
## high order finite differences
function FiniteDifference(x::Vector{Float64}, y::Vector{Float64}, d_order::Int64)::Vector{Float64}
    local dy
    #
    ## derivative parameters for forward and backwards differences
    assymetrical_first_derivative_parameters = [[-1,       1],
                                                [-1.5,     2, -0.5], 
                                                [-11/6,    3, -3/2,  1/3],
                                                [-25/12,   4, -3,	4/3,  -0.25],
                                                [-137/60,  5, -5,    10/3, -5/4,1 /5],
                                                [-49/20,   6, -15/2, 20/3, -15/4, 6/5,  -1/6],
                                                [-363/140, 7, -21/2, 35/3, -35/4, 21/5, -7/6,  1/7],
                                                [-761/280, 8, -14,   56/3, -35/2, 56/5, -14/3, 8/7, -1/8]]
    assymetrical_second_derivative_parameters = [[1,          -2,      1],	 	 	 	 	 	 
                                                 [2,          -5,      4,      -1],	 	 	 	 	 
                                                 [35/12,      -26/3,   19/2,   -14/3,    11/12], 	 	 	 
                                                 [15/4,       -77/6,   107/6,  -13,      61/12, -5/6],	 	 	 
                                                 [203/45,     -87/5,   117/4,  -254/9,   33/2,  -27/5,   137/180],	 	 
                                                 [469/90,     -223/10, 879/20, -949/18,  41,    -201/10, 1019/180, -7/10],
                                                 [29531/5040, -962/35, 621/10, -4006/45, 691/8, -282/5,  2143/90,  -206/35, 363/560]]
    assymetrical_third_derivative_parameters = [[-1,       3,      -3,         1],
                                                [-5/2,     9,      -12,        7,       -3/2],
                                                [-17/4,    71/4,   -59/2,      49/2,    -41/4,    7/4],
                                                [-49/8,    29,     -461/8,     62,      -307/8,   13,      -15/8],
                                                [-967/120, 638/15, -3929/40,   389/3,   -2545/24, 268/5,   -1849/120, 29/15],
                                                [-801/80,  349/6,  -18353/120, 2391/10, -1457/6,  4891/30, -561/8,    527/30, -469/240]]
    assymetrical_fourth_derivative_parameters =[[1,       -4,       6,        -4,      1]
                                                [3,       -14,      26,       -24,     11,       -2]
                                                [35/6,    -31,      137/2,    -242/3,  107/2,    -19,      17/6]
                                                [28/3,    -111/2,   142,      -1219/6, 176,      -185/2,   82/3,    -7/2]
                                                [1069/80, -1316/15, 15289/60, -2144/5, 10993/24, -4772/15, 2803/20, -536/15, 967/240]]
    #
    assymetrical_derivative_parameters = [assymetrical_first_derivative_parameters,
                                          assymetrical_second_derivative_parameters,
                                          assymetrical_third_derivative_parameters,
                                          assymetrical_fourth_derivative_parameters]
    #
    ## central difference derivative parameters
    central_first_derivative_parameters  = [[   -1/2,      0,   1/2],
                                            [   1/12,   -2/3,	   0,  2/3,	-1/12,],
                                            [  -1/60,   3/20,	-3/4,	 0,	  3/4, -3/20, 1/60],
                                            [  1/280, -4/105,	 1/5, -4/5,	    0,	 4/5, -1/5,	4/105, -1/280]]
    central_second_derivative_parameters = [[      1,     -2,	1],				
                                            [  -1/12,    4/3,   -5/2, 4/3,	  -1/12],			
                                            [   1/90,  -3/20,  3/2, -49/18, 3/2,	   -3/20, 1/90],		
                                            [ -1/560,  8/315, -1/5, 8/5,	  -205/72,	8/5,  -1/5, 8/315, -1/560]]
    central_third_derivative_parameters  = [[   -1/2,      1,    0,       -1,    1/2],			
                                            [    1/8,     -1,    13/8,     0,     -13/8,	1,	    -1/8],		
                                            [ -7/240,   3/10, -169/120, 61/30, 0,	    -61/30, 169/120, -3/10, 7/240]]	
    central_fourth_derivative_parameters = [[      1,     -4,	     6,	     -4,	      1],
	                                        [   -1/6,      2,	 -13/2,	   28/3,	  -13/2,	       2,	  -1/6],
	                                        [  7/240,   -2/5,	169/60,	-122/15,	   91/8,	 -122/15,	169/60,	-2/5,  7/240]]
    central_fifth_derivative_parameters  = [[   -1/2,      2,	  -5/2,	      0,	    5/2,	      -2,	   1/2],
	                                        [    1/6,   -3/2,	  13/3,	  -29/6,	      0,	    29/6,	 -13/3,   3/2,  -1/6],
	                                        [-13/288,  19/36,	-87/32,	   13/2,	-323/48,	       0,	323/48, -13/2, 87/32,	-19/36,	13/288]]
    central_sixth_derivative_parameters  = [[	   1,     -6,	    15,	    -20,	     15,	      -6,	     1],	
	                                        [   -1/4,      3,	   -13,	     29,	  -75/2,	      29,	   -13,	   3,  -1/4],	
	                                        [ 13/240, -19/24,	 87/16,	  -39/2,	  323/8,	-1023/20,	 323/8, -39/2, 87/16,	-19/24,	13/240]]
    #
    central_derivative_parameters = [central_first_derivative_parameters,
                                     central_second_derivative_parameters,
                                     central_third_derivative_parameters,
                                     central_fourth_derivative_parameters,
                                     central_fifth_derivative_parameters,
                                     central_sixth_derivative_parameters]
    #
    function forward_difference(x, f, idx, derivative_parameters)
        #
        h = mean(diff(x))
        #
        df = 0.0
        #
        for i=1:lastindex(derivative_parameters)
            df += derivative_parameters[i] * f[idx + i - 1]
        end
        #

        return df / h^(d_order)
    end
    #  
    function backward_difference(x, f, idx, derivative_parameters)
        #
        h = mean(diff(x))
        #
        df = 0.0
        #
        for i=1:lastindex(derivative_parameters)
            df += (-1)^(d_order) * derivative_parameters[i] * f[idx - i + 1]
        end        
        #
        return df / h^(d_order)
    end
    #   
    function central_difference(x,f,idx,derivative_parameters)
        #
        h = mean(diff(x))
        #
        df = 0.0
        #
        n = Int((length(derivative_parameters)-1)/2)
        #
        j = 0
        for i in -n:n
            j += 1
            df += derivative_parameters[j] * f[idx + i]
        end
        #
        return df / h^(d_order)
    end    
    #
    function max_finite_diff_orders(grid, assymetric_params, central_params)
        n = length(grid)
        #
        max_forward_orders  = zeros(Int, n)
        max_central_orders  = zeros(Int, n)
        #
        ## Forward and Backwards Differences
        for i in 1:n
            for order in 1:length(assymetric_params)
                if i + length(assymetric_params[order]) - 1 <= n
                    max_forward_orders[i] = order
                else
                    break
                end
            end
        end
        #
        ## Central differences
        for i in 1:n
            for order in 1:length(central_params)
                half_width = div(length(central_params[order]) - 1, 2)
                if (i - half_width >= 1) && (i + half_width <= n)
                    max_central_orders[i] = order
                else
                    break
                end
            end
        end
        #
        return max_forward_orders, reverse(max_forward_orders), max_central_orders
    end
    #
    ## initialise derivative vector
    dy = zeros(Float64,length(x))
    #
    ## maximum order for the requested derivative
    max_len = length(assymetrical_derivative_parameters[d_order][end])
    max_forward_order, max_backward_order, max_central_order = max_finite_diff_orders(x, assymetrical_derivative_parameters[d_order], central_derivative_parameters[d_order])
    #
    ## asymmetrical finite differences for the first and last 3 points up to maximum order
    for idx=1:3
        #
        ## maximum order for the asymmetrical finite differences
        F_order = max_forward_order[idx]
        B_order = F_order #max_backward_order[idx]
        len = Int( (length(assymetrical_derivative_parameters[d_order][F_order]) - 1))
        #
        ## forwards differences for first 3 points on grid
        dy[idx] = forward_difference(x[ idx : idx + len], y[ idx : idx + len], 1, assymetrical_derivative_parameters[d_order][F_order])
        #
        ## backwards differences for last 3 points on grid
        dy[end - idx + 1] = backward_difference(x[ end - idx + 1 - len : end - idx + 1 ], y[ end - idx + 1 - len : end - idx + 1 ], len + 1, assymetrical_derivative_parameters[d_order][B_order])
    end
    #
    function compute_central_difference(x, y, idxs)
        central_len = length(x[4 : end - 3])
        #
        if central_len % 2 == 1
            central_idx = Int((central_len - 1)/2)
        else
            central_idx = Int(central_len/2)
        end
        #
        dy = zeros(Float64,central_len)
        #
        for idx in idxs
            C_order = Int(max_central_order[idx])
            num_points = Int( (length(central_derivative_parameters[d_order][C_order]) - 1) / 2)
            #
            idx_L = 4 + idx - 1
            idx_R = lastindex(x) - 3 - idx + 1
            #
            ## first half
            dy[idx_L] = central_difference(x[idx_L - num_points: idx_L + num_points], y[idx_L - num_points: idx_L + num_points], num_points + 1, central_derivative_parameters[d_order][C_order])
            #
            ## last half
            dy[idx_R] = central_difference(x[idx_R - num_points: idx_R + num_points], y[idx_R - num_points: idx_R + num_points], num_points + 1, central_derivative_parameters[d_order][C_order])
        end
        #
        if central_len % 2 == 1
            idx_C      = start_idx + 3 + central_idx
            C_order    = Int(max_central_order[idx_C])
            num_points = Int( (length(central_derivative_parameters[d_order][C_order]) - 1) / 2)
            dy[idx_C]  = central_difference(x[idx_C - num_points: idx_C +  num_points], y[idx_C -  num_points: idx_C +  num_points],  num_points + 1,  central_derivative_parameters[d_order][C_order])
        end
        #
        return dy
    end
    #
    ## multi-threading
    # if Threads.nthreads() > 1
    #     idxs = collect(4:lastindex(x)-3)
    #     num_threads = min(Threads.nthreads(), length(idxs))
    #     chunk_size = ceil(Int, length(idxs) / num_threads)
    #     #
        # @threads for i in 1:num_threads
        #     start_idx = idxs[(i - 1) * chunk_size + 1]
        #     end_idx = Int(min(i * chunk_size, length(idxs)))
        #     end_idx = idxs[end_idx]
        #     #
        #     if start_idx <= end_idx
        #         for idx=start_idx:end_idx
        #             local C_order, num_points
        #             #
        #             C_order = Int(max_central_order[idx])
        #             num_points = Int( (length(central_derivative_parameters[d_order][C_order]) - 1) / 2)
        #             #
        #             dy[idx] = central_difference(x[idx - num_points: idx + num_points], y[idx - num_points: idx + num_points], num_points + 1, central_derivative_parameters[d_order][C_order])
        #         end
        #     end
        # end
    # else
    #
    ## central differences for middle points
    central_len = length(x[4:lastindex(x)-3])
    if central_len % 2 == 1
        central_idx = Int((central_len - 1)/2)
    else
        central_idx = Int(central_len/2)
    end
    #
    for idx=1:central_idx
        C_order = Int(max_central_order[4 + idx - 1])
        num_points = Int( (length(central_derivative_parameters[d_order][C_order]) - 1) / 2)
        #
        idx_L = 4 + idx - 1
        idx_R = lastindex(x) - 3 - idx + 1
        #
        ## first half
        dy[idx_L] = central_difference(x[idx_L - num_points: idx_L + num_points], y[idx_L - num_points: idx_L + num_points], num_points + 1, central_derivative_parameters[d_order][C_order])
        #
        ## last half
        dy[idx_R] = central_difference(x[idx_R - num_points: idx_R + num_points], y[idx_R - num_points: idx_R + num_points], num_points + 1, central_derivative_parameters[d_order][C_order])
    end

    if central_len % 2 == 1
        idx_C      = 4 + central_idx
        C_order    = Int(max_central_order[idx_C])
        num_points = Int( (length(central_derivative_parameters[d_order][C_order]) - 1) / 2)
        dy[idx_C]  = central_difference(x[idx_C - num_points: idx_C +  num_points], y[idx_C -  num_points: idx_C +  num_points],  num_points + 1,  central_derivative_parameters[d_order][C_order])
    end
    # end
    #
    return dy
end
#
function FiniteDifferenceSP(x::Vector{Float64}, y::Vector{Float64}, idx::Int64, d_order::Int64)::Float64
    local dy
    #
    ## derivative parameters for forward and backwards differences
    assymetrical_first_derivative_parameters = [[-1,       1],
                                                [-1.5,     2, -0.5], 
                                                [-11/6,    3, -3/2,  1/3],
                                                [-25/12,   4, -3,	4/3,  -0.25],
                                                [-137/60,  5, -5,    10/3, -5/4,1 /5],
                                                [-49/20,   6, -15/2, 20/3, -15/4, 6/5,  -1/6],
                                                [-363/140, 7, -21/2, 35/3, -35/4, 21/5, -7/6,  1/7],
                                                [-761/280, 8, -14,   56/3, -35/2, 56/5, -14/3, 8/7, -1/8]]
    assymetrical_second_derivative_parameters = [[1,          -2,      1],	 	 	 	 	 	 
                                                 [2,          -5,      4,      -1],	 	 	 	 	 
                                                 [35/12,      -26/3,   19/2,   -14/3,    11/12], 	 	 	 
                                                 [15/4,       -77/6,   107/6,  -13,      61/12, -5/6],	 	 	 
                                                 [203/45,     -87/5,   117/4,  -254/9,   33/2,  -27/5,   137/180],	 	 
                                                 [469/90,     -223/10, 879/20, -949/18,  41,    -201/10, 1019/180, -7/10],
                                                 [29531/5040, -962/35, 621/10, -4006/45, 691/8, -282/5,  2143/90,  -206/35, 363/560]]
    assymetrical_third_derivative_parameters = [[-1,       3,      -3,         1],
                                                [-5/2,     9,      -12,        7,       -3/2],
                                                [-17/4,    71/4,   -59/2,      49/2,    -41/4,    7/4],
                                                [-49/8,    29,     -461/8,     62,      -307/8,   13,      -15/8],
                                                [-967/120, 638/15, -3929/40,   389/3,   -2545/24, 268/5,   -1849/120, 29/15],
                                                [-801/80,  349/6,  -18353/120, 2391/10, -1457/6,  4891/30, -561/8,    527/30, -469/240]]
    assymetrical_fourth_derivative_parameters =[[1,       -4,       6,        -4,      1]
                                                [3,       -14,      26,       -24,     11,       -2]
                                                [35/6,    -31,      137/2,    -242/3,  107/2,    -19,      17/6]
                                                [28/3,    -111/2,   142,      -1219/6, 176,      -185/2,   82/3,    -7/2]
                                                [1069/80, -1316/15, 15289/60, -2144/5, 10993/24, -4772/15, 2803/20, -536/15, 967/240]]
    #
    assymetrical_derivative_parameters = [assymetrical_first_derivative_parameters,
                                          assymetrical_second_derivative_parameters,
                                          assymetrical_third_derivative_parameters,
                                          assymetrical_fourth_derivative_parameters]
    #
    ## central difference derivative parameters
    central_first_derivative_parameters  = [[   -1/2,      0,   1/2],
                                            [   1/12,   -2/3,	   0,  2/3,	-1/12,],
                                            [  -1/60,   3/20,	-3/4,	 0,	  3/4, -3/20, 1/60],
                                            [  1/280, -4/105,	 1/5, -4/5,	    0,	 4/5, -1/5,	4/105, -1/280]]
    central_second_derivative_parameters = [[      1,     -2,	1],				
                                            [  -1/12,    4/3,   -5/2, 4/3,	  -1/12],			
                                            [   1/90,  -3/20,  3/2, -49/18, 3/2,	   -3/20, 1/90],		
                                            [ -1/560,  8/315, -1/5, 8/5,	  -205/72,	8/5,  -1/5, 8/315, -1/560]]
    central_third_derivative_parameters  = [[   -1/2,      1,    0,       -1,    1/2],			
                                            [    1/8,     -1,    13/8,     0,     -13/8,	1,	    -1/8],		
                                            [ -7/240,   3/10, -169/120, 61/30, 0,	    -61/30, 169/120, -3/10, 7/240]]	
    central_fourth_derivative_parameters = [[      1,     -4,	     6,	     -4,	      1],
	                                        [   -1/6,      2,	 -13/2,	   28/3,	  -13/2,	       2,	  -1/6],
	                                        [  7/240,   -2/5,	169/60,	-122/15,	   91/8,	 -122/15,	169/60,	-2/5,  7/240]]
    central_fifth_derivative_parameters  = [[   -1/2,      2,	  -5/2,	      0,	    5/2,	      -2,	   1/2],
	                                        [    1/6,   -3/2,	  13/3,	  -29/6,	      0,	    29/6,	 -13/3,   3/2,  -1/6],
	                                        [-13/288,  19/36,	-87/32,	   13/2,	-323/48,	       0,	323/48, -13/2, 87/32,	-19/36,	13/288]]
    central_sixth_derivative_parameters  = [[	   1,     -6,	    15,	    -20,	     15,	      -6,	     1],	
	                                        [   -1/4,      3,	   -13,	     29,	  -75/2,	      29,	   -13,	   3,  -1/4],	
	                                        [ 13/240, -19/24,	 87/16,	  -39/2,	  323/8,	-1023/20,	 323/8, -39/2, 87/16,	-19/24,	13/240]]
    #
    central_derivative_parameters = [central_first_derivative_parameters,
                                     central_second_derivative_parameters,
                                     central_third_derivative_parameters,
                                     central_fourth_derivative_parameters,
                                     central_fifth_derivative_parameters,
                                     central_sixth_derivative_parameters]
    #
    #
    # Helper functions for the three types of finite differences (same as original code)
    function forward_difference(f, derivative_parameters, h)
        df = 0.0
        for i=1:lastindex(derivative_parameters)
            df += derivative_parameters[i] * f[i]
        end
        return df / h^(d_order)
    end
    function backward_difference(f, derivative_parameters, h)
        df = 0.0
        for i=1:lastindex(derivative_parameters)
            df += (-1)^(d_order) * derivative_parameters[i] * f[end - i + 1]
        end        
        return df / h^(d_order)
    end
    function central_difference(f, derivative_parameters, h)
        df = 0.0
        n = Int((length(derivative_parameters)-1)/2)
        j = 0
        for i in -n:n
            j += 1
            df += derivative_parameters[j] * f[n + 1 + i]
        end
        return df / h^(d_order)
    end    
    #
    ##
    #
    # Main logic to compute the derivative at a single point
    h = mean(diff(x))
    n = length(x)
    
    # Check for valid d_order
    if d_order < 1 || d_order > length(assymetrical_derivative_parameters)
        throw(ArgumentError("Derivative order not supported."))
    end
    
    # Try Central Difference first as it is the most accurate
    # Find the best order for central difference at this point
    central_params = central_derivative_parameters[d_order]
    best_central_order = 0
    for order = 1:length(central_params)
        half_width = div(length(central_params[order]) - 1, 2)
        if (idx - half_width >= 1) && (idx + half_width <= n)
            best_central_order = order
        else
            break
        end
    end

    if best_central_order > 0

        # Compute central difference
        params = central_params[best_central_order]
        half_width = div(length(params) - 1, 2)
        y_slice = y[idx - half_width : idx + half_width]
        return central_difference(y_slice, params, h)
    else
        # Fall back to forward or backward difference
        assym_params = assymetrical_derivative_parameters[d_order]
        
        # Check for forward difference
        best_forward_order = 0
        for order = 1:length(assym_params)
            num_points = length(assym_params[order])
            if idx + num_points - 1 <= n
                best_forward_order = order
            else
                break
            end
        end

        if best_forward_order > 0
            # Compute forward difference
            params = assym_params[best_forward_order]
            num_points = length(params)
            y_slice = y[idx : idx + num_points - 1]
            return forward_difference(y_slice, params, h)
        else
            # Must use backward difference
            best_backward_order = 0
            for order = 1:length(assym_params)
                num_points = length(assym_params[order])
                if idx - num_points + 1 >= 1
                    best_backward_order = order
                else
                    break
                end
            end
            
            if best_backward_order > 0
                # Compute backward difference
                params = assym_params[best_backward_order]
                num_points = length(params)
                y_slice = y[idx - num_points + 1 : idx]
                return backward_difference(y_slice, params, h)
            else
                throw(ArgumentError("Cannot compute derivative at this point with the available grid points."))
            end
        end
    end
end
#
## cumulative integral
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
## simpson integration
@views function simps(f::Function, x::AbstractRange) #A range with odd number
    h = step(x)
    I= h/3*(f(x[1])+2*sum(f,x[3:2:end-2])+4*sum(f,x[2:2:end-1])+f(x[end]))
    return I
end
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

#
function is_triangular(x::Int)
    if x < 0
        return false
    end
    
    # Calculate 1 + 8x
    val = 1 + 8x

    # Check if val is a perfect square
    s = isqrt(val)
    if s * s != val
        return false
    end

    # Check if the result of N is an integer
    return (-1 + s) % 2 == 0
end
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GRID FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
function gridSplitter(Ninit,Nsplit)
    N = Ninit
    for i=1:Nsplit
        N = 2*N -1
    end
    #
    return N
end
#
function extend_grid_left(grid, to_value)
    # Determine the current step size
    step = grid[2] - grid[1]
    #
    new_grid = []
    current_value = grid[1] - step
    #
    while current_value > to_value
        push!(new_grid, current_value)
        current_value -= step
    end
    push!(new_grid,to_value)
    return vcat(reverse(new_grid), grid)
end
#
function smoothgrid(rmin,rmax,rgrid,nsplit,nsL,nsR)
    #
    ## extend the rgrid to zero with the same cadence
    rgrid = extend_grid_left(rgrid, 0.0)
    r0 = LinRange(rmin, rgrid[1]-(1000+rgrid[1])/nsL, nsL)
    r1 = LinRange(rgrid[1], rgrid[end], gridSplitter(size(rgrid)[1],nsplit))
    r2 = LinRange(rgrid[end]+(1000-rgrid[end])/nsR,rmax, nsR)
    rsolve = []
    append!(rsolve,r0)
    append!(rsolve,r1)
    append!(rsolve,r2)
    rsolve = unique(sort(rsolve))
    rsolve = map(Float64, rsolve)
    #
    ∆0 = r0[2] - r0[1]
    ∆1 = r1[2] - r1[1]
    ∆2 = r2[2] - r2[1]
    # println(∆0," ",∆1," ",∆2)
    #
    gamma = 0.02 #-log(1e-8)/(rsolve[flagidx-1]-rsolve[zeroidx])
    s1 = sigmoid.(rsolve,gamma,(rmin + rgrid[1])/2)
    s2 = sigmoid.(rsolve,gamma,(rmax + rgrid[end])/2)
    ∆ = []
    for i=1:lastindex(rsolve)
        ∆a  =  s1[i]*∆1 + (1 - s1[i])*∆0
        ∆i  =  s2[i]*∆2 + (1-s2[i])*∆a  
        push!(∆,∆i)
    end
    #
    diffs = Spline1D(rsolve,∆)
    #
    rL = r1[1] - sum(∆[1:10001])
    #
    r = zeros(Float64, length(rsolve))
    r[1] = rL
    for idx=2:length(rsolve)
        r[idx] = r[idx-1] + ∆[idx]
    end
    return r, ∆
end
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PROPERTY FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
function sigmoid(x,g,x0)
    return (1 + exp(-g*(x-x0)))^(-1)
end
#
function double_sigmoid(x,y,g,x0)
    #
    P = []
    z = zeros(Float64, length(x))

    s1 = sigmoid.(x,g[1],x0[1])
    s2 = sigmoid.(x,g[2],x0[2])

    a  = s1.*y .- (s1 .- 1).*z
    b  = s2.*z .- (s2 .- 1).*a

    return b
end
#
function repulsive(r,params...)
    coefficients = params
    #
    V = 0
    #
    for i=1:lastindex(coefficients) 
        V += coefficients[i]/r^(i-1)
    end
    #
    return V
end
#
function morse_oscillator(r,Ve,re,a,Ae) 
    exponent = a
    return Ve+(Ae-Ve)*(1 - exp(-exponent*(r-re)))^(2)
end
#
function EMO(R, VE, RE, AE, RREF, PL, PR, NL, NR, params...)    
    param = params
    #
    exponent = 0
    #
    if RREF == -1.0; RREF = RE; end
    #
    if R<RE
        y = (R^PL-RREF^PL)/(R^PL+RREF^PL) #Surkus Variable
        exponent = 0
        for i=1:round(Int,NL)+1
            exponent += y^(i-1)*param[i]
        end
        return VE+(AE-VE)*(1 - exp(-exponent*(R-RE)))^(2)
    end
    if R>=RE
        y = (R^PR-RREF^PR)/(R^PR+RREF^PR) #Surkus Variable
        exponent = 0
        for i=1:round(Int,NR)+1
            exponent += y^(i-1)*param[i]
        end
        return VE+(AE-VE)*(1 - exp(-exponent*(R-RE)))^(2)
    end
end
#
function surkus_expansion(R,RE,RREF,p,expansionParams)
    if RREF==-1.0; RREF=RE; end
    y = (R^(p)-RREF^(p))/(R^(p)+RREF^(p))
    expansion = 0
    for (i, param) in enumerate(expansionParams)
        expansion += param*y^(i-1)
    end
    return expansion
end
#
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
function DC_gaussian(r, s, r0, rref, p, expansionParams)
    f = surkus_expansion(r,r0,rref,p,expansionParams)
    return f*exp(-(r-r0)^2/(2s^2))/(s*sqrt(2*pi))
end
#
function DC_lorentzian(r, w, r0, rref, p, expansionParams)          # Lorentzian 
    f = surkus_expansion(r,r0,rref,p,expansionParams)
    return 0.5*f*(w/((r - r0)^2 .+ w^2))
end
#
function DC_laplacian(r, b, r0, rref, p, expansionParams)            # Laplacian
    f = surkus_expansion(r,r0,rref,p,expansionParams)
	return f*(pi/(4*b))*exp(-( (abs(r-r0))/b ))
end
#
function mixAng_lo(r,w,rc)                     # Lorentizan mixing angle
    return (pi/4)+(atan((r-rc)/w)/2)
end
#
function mixAng_la(r,b,rc)                       # Laplacian mixing angle
	if r<rc
		return (pi/4)*exp((r-rc)/b)
	else
		return (pi/2)-(pi/4)*exp(-(r-rc)/b)
	end
end
#
function ref_sigmoid(x, ref, a, x0)
    return  2*ref/(1+exp(a*(x-x0)))
end
#
function mix(r, g0, r0, amp, a, m, sub_types)
    #
    ## evalute functional symbols
    f1 = eval(Symbol(sub_types[1]))
    f2 = eval(Symbol(sub_types[2]))
    #
    ## compute gamma
    if a != 0.0
        g = ref_sigmoid(r,g0,a,r0)
    else
        g = g0
    end
    #
    ## compute mixture
    return m * f2(r, g, r0, amp) + (1 - m) * f1(r, g, r0, amp)
end
#
# function mix_perturbed(r, g0, r0, amp, a, m, ptrb_and_subtype_params...)
#     sub_types = ptrb_and_subtype_params[end]
#     p_ptrb = collect(ptrb_and_subtype_params[1:end-1])
#     #
#     ## evalute functional symbols
#     f1 = eval(Symbol(sub_types[1]))
#     f2 = eval(Symbol(sub_types[2]))
#     #
#     ## compute gamma
#     if a != 0.0
#         g = ref_sigmoid(r,g0,a,r0)
#     else
#         g = g0
#     end
#     #
#     ## compute perturbation
#     fP = eval(Symbol(sub_types[3]))
#     #
#     ## compute how many paerturbations
#     N_ptrb = Int(length(p_ptrb)/3)
#     p_ptrb = reshape(p_ptrb,N_ptrb,:)
#     #
#     P = 0.0
#     for idx=1:N_ptrb
#         gi, r0i, ampi = p_ptrb[1,idx], p_ptrb[2,idx], p_ptrb[3,idx]
#         P+=fP(r, gi, r0i, ampi)
#     end
#     #
#     ## compute mixture + perturbation
#     return abs(m) * f2(r, g, r0, amp) + (1 - abs(m)) * f1(r, g, r0, amp) + P
# end
function mix_perturbed(r, g0, r0, amp, a, m, ptrb_and_subtype_params...)
    sub_types = ptrb_and_subtype_params[end]
    p_ptrb = collect(ptrb_and_subtype_params[1:end-1])
    #
    ## evalute functional symbols
    f1 = eval(Symbol(sub_types[1]))
    f2 = eval(Symbol(sub_types[2]))
    #
    ## compute gamma
    if a != 0.0
        g = ref_sigmoid(r,g0,a,r0)
    else
        g = g0
    end
    #
    ## compute perturbation
    fP = eval(Symbol(sub_types[3]))
    #
    ## compute how many paerturbations
    N_ptrb = Int(length(p_ptrb)/3)
    p_ptrb = reshape(p_ptrb,:,N_ptrb)
    #
    P = 0.0
    if N_ptrb == 1
        gi, r0i, ampi = p_ptrb[1], p_ptrb[2], p_ptrb[3]
        P+=fP(r, gi, r0i, ampi)
    else
        for idx=1:N_ptrb
            gi, r0i, ampi = p_ptrb[1,idx], p_ptrb[2,idx], p_ptrb[3,idx]
            P+=fP(r, gi, r0i, ampi)
        end
    end
    #
    ## compute mixture + perturbation
    return abs(m) * f2(r, g, r0, amp) + (1 - abs(m)) * f1(r, g, r0, amp) + P
end
#
function lorentzian(r,w,rc,amp)
    return amp*0.5*(w/((r - rc)^2 + w^2))
end
#
function duo_lorentzian(r,w,rc,a0)
    return 2*(a0/pi)*(w/(4*(r-rc)^2 + w^2))
end
#
function laplace(r,b,rc)                                             # Laplacian
	return (pi/(4*b))exp(-( (abs.(r-rc))/b ))
end
#
function gaussian(r,w,r0,amp)
    gauss= (0.5*amp/w)*exp(-log(2)*((r-r0)/w)^2)
    if isnan(gauss)
        return 0.0
    else
        return gauss
    end
end 
#
function lor_lap(r,w,rc)
    #
    ## compute geom average lorentzian laplacian mixing angle
    beta_geomavg = mixAng_lor_lap.(r,w,rc)
    # println(beta_geomavg)
    #
    ## compute the derivative of the mixing angle to obtain the NAC
    return sign(w)*FiniteDifference(collect(r),beta_geomavg,1)
end
#
function mixAng_duo_lorentzian(r,w,rc,a0)                             # Lorentizan mixing angle
	return (amp*pi/4)+(amp*atan((r-rc)/w)/2)
end
#
function mixAng_lorentzian(r,w,rc,amp)                             # Lorentizan mixing angle
	return (amp*pi/4)+(amp*atan((r-rc)/w)/2)
end
#
function mixAng_laplacian(r,b,rc)                              # Laplacian mixing angle
	if r<rc
		return (pi/4)*exp((r-rc)/b)
	else
		return (pi/2)-(pi/4)*exp(-(r-rc)/b)
	end
end
#
function mixAng_gaussian(xi,x0,A,µ,d)
        return (A*d*sqrt(pi)/sqrt(2))*(erf( (xi-µ)/(d*sqrt(2)) ) + erf( (x0+µ)/(d*sqrt(2)) ))
end
#
# See DOI: 10.1039/D2CP03051A for more details
# function mixAng_avg(r,w,rc) # geometrical averaged Lorent-Laplace mixing angle
#     w = abs(w)      # absolute Lorentz HWHM to stop bound errors in optimization
#     b = 1.397*w     # relation of w to Laplace parameter through maximal overlap 
# 	if r<rc	         
# 		return          (1/2)*asin( sqrt(  sin(2*mixAng_lo(r,w,rc))  *  sin(2*mixAng_la(r,b,rc))  ) )
# 	else
# 		return (pi/2) - (1/2)*asin( sqrt(  sin(2*mixAng_lo(r,w,rc))  *  sin(2*mixAng_la(r,b,rc))  ) )
# 	end
# end
#
function mixAng_lor_lap(r,w,rc)   # geometrical averaged Lorent-Laplace mixing angle
    w = abs(w)      # absolute Lorentz HWHM to stop bound errors in optimization
    b = 1.397*w     # relation of w to Laplace parameter through maximal overlap 
	if r<rc	         
		return          (1/2)*asin( sqrt(  sin(2*mixAng_lo(r,w,rc))  *  sin(2*mixAng_la(r,b,rc))  ) )
	else
		return (pi/2) - (1/2)*asin( sqrt(  sin(2*mixAng_lo(r,w,rc))  *  sin(2*mixAng_la(r,b,rc))  ) )
	end
end
#
function DC_beta(R, V1, V2, BETA, ARGS...) #GAMMA, RC
    mixAngle = map(x -> eval(Symbol("mixAng_"*BETA))(x,ARGS...), R)
    PotenDiffOver2 = (V2.-V1)./2
	return PotenDiffOver2.*tan.(mixAngle.*2)
end
#
function discriminant(a,b,c,d)
    lp = ((a.+d)./2).+0.5*sqrt.((a.-d).^2 .+(b.*c))
    lm = ((a.+d)./2).-0.5*sqrt.((a.-d).^2 .+(b.*c)) 
    return lm, lp
end
#
function discriminant_real_symmetric(a,b,d)
    lp = ((a.+d)./2).+0.5*sqrt.((a.-d).^2 .+ (4 .* b))
    lm = ((a.+d)./2).-0.5*sqrt.((a.-d).^2 .+ (4 .* b)) 
    return lm, lp
end
#
function COUPLED_PEC(r, LVAL, RVAL, TYPES, Np, DIMENSION; betaCouple=false, crossing = false)
    #
    ## index mapping function
    function index_map(i,j,N)
        if i == j
            return i
        else
            return N + (i-1)*N - ((i-1)*i/2) + (j-i)
        end
    end
    #
    function generate_ij_pairs(DIMENSION::Int)
        #
        ## Part 1: Diagonal elements (i,i)
        diagonal_pairs = [(i, i) for i = 1:DIMENSION]
        #
        ## Part 2: Off-diagonal elements (i,j) where i < j, ordered by i then j
        off_diagonal_pairs = Tuple{Int,Int}[] # Initialize an empty array of tuples
        for i = 1:DIMENSION
            for j = (i+1):DIMENSION # j starts from i+1 to ensure i < j
                push!(off_diagonal_pairs, (i, j))
            end
        end
        #
        ## Combine the two parts. The `vcat` function concatenates arrays
        return vcat(diagonal_pairs, off_diagonal_pairs)
    end
    #
    RC = -1.0
    #
    ## find the component the user wishes to compute
    COMPONENT = Int(RVAL[end])
    #
    ## slice the parameter arrays according to the number of parameters in Np
    L = []
    R = []
    i = 1
    for N in Np
        f = Int(i + N - 1)
        #
        ## populate array slices
        push!(L, LVAL[i:f])
        push!(R, RVAL[i:f])
        #
        ## update start index
        i = f + 1
    end
    V = zeros(Float64, length(r), DIMENSION, DIMENSION)
    #
    mu_idx = generate_ij_pairs(DIMENSION)
    #
    for (mu, pair) in enumerate(mu_idx)
        #
        i, j = pair
        #
        if (i!=j)&(betaCouple==true)
            V1d = V[:,i,i]
            V2d = V[:,j,j]
            #
            ## find crossing point
            if crossing == false
                spl = Spline1D(r, map(i->V2d[i]-V1d[i],1:length(r)))
                root = Dierckx.roots(spl)
                if !isempty(root)
                    RC = root[1]
                else
                    RC = R[mu][2]
                end
            else
                RC = crossing
            end
            #
            R[mu][2] = RC
            #
            y = DC_beta(r,
                        V1d,
                        V2d,
                        TYPES[mu],
                        R[mu]...)
        else
            _, y = ComputeProperty_viaParameters(r, 
                                                TYPES[mu],
                                                L[mu], 
                                                R[mu],
                                                "pec",
                                                ["angstrom","cm-1"],
                                                ("N/A","N/A","N/A"),
                                                1)
        end
        #
        V[:,i,j] .= y
    end
    #
    ## diagonalise the potential matrix
    if DIMENSION == 2  
        V1a, V2a = discriminant_real_symmetric(V[:,1,1],V[:,1,2],V[:,2,2])
        #
        adi = [V1a, V2a]
        #
        return adi[COMPONENT], RC
    else
        Va = map(x -> eigen(V[x,:,:]).values[COMPONENT], collect(1:length(r)))
        #
        return Va, RC
    end
end

#
# function sigmoid_polynom_decay_24(r, gamma0, r0, )
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ METHOD FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
function ComputeProperty(self; custom_grid = false, evolution_grid = false)
    #
    ## define grid interpolated bond lengths
    if custom_grid == false
        r  = LinRange(Calculation["grid"].range[1],
                      Calculation["grid"].range[2],
                      Calculation["grid"].npoints)
    else
        r = custom_grid
    end
    #
    if factor == 0.0
        return r, zeros(Float64,length(r))
    end
    #
    if self.type == "grid"                                      # grid based PEC
        local r, spline, f_interpolated
        #
        # DEFAULT IS CUBIC FOR NOW *********************************************
		## spline the property onto the defined grid
        # if occursin("linear",Calculation["grid"].interpolation_type)
        #     spline = LinearInterpolation(self.Lval, self.Rval)  
        # elseif occursin("quadratic",Calculation["grid"].interpolation_type)
        #     spline = interpolate((self.Lval,), self.Rval, Gridded(Quadratic()))
        # elseif occursin("cubic",Calculation["grid"].interpolation_type)
        #     spline = Spline1D(self.Lval, self.Rval)                
            #spline = CubicSplineInterpolation(self.Lval, self.Rval)
        # elseif occursin("quintic",Calculation["grid"].interpolation_type)
        #     spline = Fun(self.Lval, self.Rval, kind="quintic")
        # end              
        #
        ## convert units
        x, f = unitConversion(self.Lval, self.Rval, self.obj_type, self.units)
		#
        ## compute the interpolated object curve on the bonds grid
		spline = Spline1D(x, f, k=3, bc="extrapolate")
        f_interpolated = spline(r)
        #
        ## if evolution is being done then the interpolated grid needs to go to 
        ## 0 either side of the defined grid.
        if evolution_grid == true
            f_interpolated = double_sigmoid(r,f_interpolated,[100,100],[Calculation["grid"].range[1],Calculation["grid"].range[2]])
        end
        # #
        # ## convert units
        # r, f_interpolated = unitConversion(r, f_interpolated, self.obj_type, self.units)
        return r, f_interpolated.*self.factor
    #
    elseif any(occursin.(["coupled_pec","coupled-pec"], lowercase(self.type)))
        #
        coupled_beta = false
        #
        if occursin("beta",lowercase(self.type)); coupled_beta = true; end
        #
        Nparameters = []
        subtypes = []
        #
        for s in self.sub_type
            #
            subtype, n = split(s, "-")
            #
            push!(subtypes, subtype)
            #
            push!(Nparameters,  parse(Float64,n))
        end
        #
        ## compute the dimension of the electronic matrix to diagonalise
        Num_of_el = length(self.sub_type)
        #
        if !is_triangular(Num_of_el)
            object = self.obj_type
            object_id = self.ID
            throw(DomainError(x, 
                "The number of coupled objects for $object $object_id ($Num_of_el) is incorrect." *
                "It must be equal to N(N+1)/2, where N is the number of coupled electronic states."
            ))
        end
        #
        DIMENSION = Int((-1+sqrt(1+8*Num_of_el))/2)
        #
        f, RC = COUPLED_PEC(collect(r), self.Lval, self.Rval, subtypes, Nparameters, DIMENSION, betaCouple = coupled_beta)
        #
        #
        self.crossing = RC
        #
        ## convert units
        r, f = unitConversion(r, f, self.obj_type, self.units)
        return r, f .* self.factor
    #
    else
        #
        ## 
        func = eval(Symbol(self.type))
        #
        ## compute curve
        if self.sub_type[1] == "mixang"
            f = func(r,self.Rval...)
        elseif all(self.sub_type .== "N/A")
            if self.type == "lor_lap"
                f = func(collect(r),self.Rval...)
            else
                f = map(x -> func(x,self.Rval...),r)
            end
        else
            f = map(x -> func(x,self.Rval...,self.sub_type),r)
        end
        #
        ## convert units
        r, f = unitConversion(r, f, self.obj_type, self.units)
        return r, f .* self.factor
    end
end
#
function ComputeProperty_viaParameters(X, ftype, Lval, Rval, obj_type, units, sub_type, factor; state = false)
    r = X
    #
    if factor == 0.0
        return r, zeros(Float64,length(r))
    end
    #
    if ftype == "grid"                                    
        local r, spline, f_interpolated
        #
        # DEFAULT IS CUBIC FOR NOW *********************************************
		## spline the property onto the defined grid
        # if occursin("linear",Calculation["grid"].interpolation_type)
        #     spline = LinearInterpolation(self.Lval, self.Rval)  
        # elseif occursin("quadratic",Calculation["grid"].interpolation_type)
        #     spline = interpolate((self.Lval,), self.Rval, Gridded(Quadratic()))
        # elseif occursin("cubic",Calculation["grid"].interpolation_type)
        #     spline = Spline1D(self.Lval, self.Rval)                
            #spline = CubicSplineInterpolation(self.Lval, self.Rval)
        # elseif occursin("quintic",Calculation["grid"].interpolation_type)
        #     spline = Fun(self.Lval, self.Rval, kind="quintic")
        # end              
        #
        ## convert units
        x, f = unitConversion(Lval, Rval, obj_type, units)
		#
        ## compute the interpolated object curve on the bonds grid
		spline = Spline1D(x, f, k=3, bc="extrapolate") 
        f_interpolated = spline(r)
        #
        return r, f_interpolated.*factor
    elseif any(occursin.(["coupled_pec","coupled-pec"], lowercase(ftype)))
        #
        coupled_beta = false
        #
        if occursin("beta",lowercase(ftype)); coupled_beta = true; end
        #
        Nparameters = []
        subtypes = []
        #
        for s in sub_type
            #
            subtype, n = split(s, "-")
            #
            push!(subtypes, subtype)
            #
            push!(Nparameters,  parse(Float64,n))
        end
        #
        ## compute the dimension of the electronic matrix to diagonalise
        Num_of_el = length(sub_type)
        #
        if !is_triangular(Num_of_el)
            object = obj_type
            object_id = ID
            throw(DomainError(x, 
                "The number of coupled objects for $object $object_id ($Num_of_el) is incorrect." *
                "It must be equal to N(N+1)/2, where N is the number of coupled electronic states."
            ))
        end
        #
        DIMENSION = Int((-1+sqrt(1+8*Num_of_el))/2)
        #
        f, RC = COUPLED_PEC(collect(r), Lval, Rval, subtypes, Nparameters, DIMENSION, betaCouple = coupled_beta)
        #
        if state != false
            Potential[state].crossing = RC
        end
        #
        ## convert units
        r, f = unitConversion(r, f, obj_type, units)
        return r, f .* factor
    else
        #
        ## 
        func = eval(Symbol(ftype))
        #
        ## compute curve
        if sub_type[1] == "mixang"
            f = func(r,Rval...)
        elseif all(sub_type .== "N/A")
            f = map(x -> func(x, Rval...),r)
        else
            f = map(x -> func(x, Rval..., sub_type), r)
        end
        #
        ## convert units
        r, f = unitConversion(r, f, obj_type, units)
        return r, f.*factor
    end
end
#
function ComputeProperty_viaParameters_SP(X, ftype, Lval, Rval, obj_type, units, sub_type, factor; state = false)
    r = X
    #
    if factor == 0.0
        return 0.0
    end
    #
    if ftype == "grid"                                    
        local r, spline, f_interpolated
        #
        # DEFAULT IS CUBIC FOR NOW *********************************************
		## spline the property onto the defined grid
        # if occursin("linear",Calculation["grid"].interpolation_type)
        #     spline = LinearInterpolation(self.Lval, self.Rval)  
        # elseif occursin("quadratic",Calculation["grid"].interpolation_type)
        #     spline = interpolate((self.Lval,), self.Rval, Gridded(Quadratic()))
        # elseif occursin("cubic",Calculation["grid"].interpolation_type)
        #     spline = Spline1D(self.Lval, self.Rval)                
            #spline = CubicSplineInterpolation(self.Lval, self.Rval)
        # elseif occursin("quintic",Calculation["grid"].interpolation_type)
        #     spline = Fun(self.Lval, self.Rval, kind="quintic")
        # end              
        #
        ## convert units
        x, f = unitConversion(Lval, Rval, obj_type, units)
		#
        ## compute the interpolated object curve on the bonds grid
		spline = Spline1D(x, f, k=3, bc="extrapolate")
        f_interpolated = spline(r)
        #
        return f_interpolated[1] * factor
    elseif any(occursin.(["coupled_pec","coupled-pec"], lowercase(ftype)))
        #
        coupled_beta = false
        #
        if occursin("beta",lowercase(ftype)); coupled_beta = true; end
        #
        Nparameters = []
        subtypes = []
        #
        for s in sub_type
            #
            subtype, n = split(s, "-")
            #
            push!(subtypes, subtype)
            #
            push!(Nparameters,  parse(Float64,n))
        end
        #
        ## compute the dimension of the electronic matrix to diagonalise
        Num_of_el = length(sub_type)
        #
        if !is_triangular(Num_of_el)
            object = obj_type
            object_id = ID
            throw(DomainError(x, 
                "The number of coupled objects for $object $object_id ($Num_of_el) is incorrect." *
                "It must be equal to N(N+1)/2, where N is the number of coupled electronic states."
            ))
        end
        #
        DIMENSION = Int((-1+sqrt(1+8*Num_of_el))/2)
        #
        f,RC = COUPLED_PEC(r, Lval, Rval, subtypes, Nparameters, DIMENSION, betaCouple = coupled_beta, crossing = Potential[state].crossing)
        #
        ## convert units
        r, f = unitConversion(r, f, obj_type, units)
        return f[1] * factor 
    else
        #
        ## 
        func = eval(Symbol(ftype))
        #
        ## compute curve
        if sub_type[1] == "mixang"
            f = func(r,Rval...)
        elseif all(sub_type .== "N/A")
            f = map(x -> func(x, Rval...),r)
        else
            f = map(x -> func(x, Rval..., sub_type), r)
        end
        #
        ## convert units
        r, f = unitConversion(r, f, obj_type, units)
        return f[1] * factor
    end
end
#
function unitConversion(x, y, obj, units)
    #
    ## define an atomic unit standard for different objects
    unitStd = Dict()
    unitStd["pec"] = "eh"
    unitStd["poten"] = "eh"
    unitStd["soc"] = "cm-1"
    unitStd["dm"]  = "ea0"
    unitStd["nac"] = "ang-1"
    unitStd["lx"] = "cm-1"
    #
    ## define unit dictionary
    unitDict = Dict()
    unitDict["ang-1"]    =                       1
    unitDict["cm-1"]     =                       1
    unitDict["angstrom"] =                       1
    unitDict["debye"]    =                       1
    unitDict["bohr"]     =          0.529177210920
    unitDict["eh"]       =       219474.6313708000
    unitDict["ea0"]      =          2.541746363812
    unitDict["erg"]      =      5.034117008194E+15
    unitDict["ev"]       =         8065.5442959967
    unitDict["kcal/mol"] =          349.7550878997
    unitDict["kj/mol"]   =           83.5934722514
    unitDict["thz"]      =           33.3564095198
    unitDict["au"]       =  unitDict[unitStd[lowercase(obj)]]
    #
    ## do the unit conversion
    x = x.*unitDict[lowercase(units[1])]
    y = y.*unitDict[lowercase(units[2])]
    #
    return x, y
end
#
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
function closest_value_index(arr, value)
    index = argmin(abs.(arr .- value))
    return index
end
#
function FiniteDiff_MatDerivative(x, M, dim, d_order)
    dM = zeros(Float64,length(x),dim,dim)
    #
    ## extract elements of U, spline them, form splined matrix and derivative
    for i=1:dim
        for j=1:dim
            #
            Mij = M[:,i,j]
            #
            dMij = FiniteDifference(x,Mij,d_order)
            #
            dM[:,i,j] = dMij
        end
    end
    #
    return dM
end
#
function SplineMat(x::AbstractVector{Float64},A::AbstractArray{Float64},r::AbstractVector{Float64})::AbstractArray{Float64}
    dim = size(A[1,:,:])[1]
    B = zeros(Float64,length(r),dim,dim)
    #
    for i=1:dim
        for j=1:dim
            splij = Spline1D(x,A[:,i,j])
            B[:,i,j] = splij(r)
        end
    end
    #
    return B
end
#
function frobeniusNorm(A)
    sqr = []
    dim = size(A)[1]
    #
    for i=1:dim
        for j=1:dim
            push!(sqr,abs(A[i,j])^2)
        end
    end
    #
    return sqrt(sum(sqr))
end
#
function Diabatic_Property_from_wavefunctions(bra, ket, adiMat)
    #
    ##
    Pd_ij = map(i -> bra[i]' * adiMat[i,:,:] * ket[i], collect(1:length(bra)))
    #
    return Pd_ij
end
#
function fit_abinitio()
    function get_loss_general(x_ai, y_ai, Lval, units, obj_type, sub_type, factor, f, fixParams, p_init, p_bounds, p)
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
        ## compute curve with new parameters
        x, y = ComputeProperty_viaParameters(x_ai, f, Lval, p, obj_type, units, sub_type, factor)
        #
        ## compute RMSE
        diff = y .- y_ai
        RMSE = sqrt(sum(diff.^2)/length(diff))
        #
        return round(RMSE, sigdigits=6)
    end
    #
    options = Optim.Options(show_trace = true)
    #
    fit_flag = false
    #
    for key in keys(abinitio)
        #
        if key[1] == "poten"
            if key[2] in Calculation["method"].states
                fit_flag = true
                #
                i = floor(Int64, key[2][1])
                j = floor(Int64, key[2][1])
            end
        else
            bra_state = key[2][1]
            ket_state = key[2][2]
            #
            if (bra_state in Calculation["method"].states)&(ket_state in Calculation["method"].states)
                fit_flag = true
                #
                i = floor(Int64, key[2][1])
                j = floor(Int64, key[2][2])
            end
        end
        #
        if fit_flag == true
            fitting_object = Hamiltonian[key]
            #
            x_ai = abinitio[key].Lval
            y_ai = abinitio[key].Rval * abinitio[key].factor
            #
            ## make mask for fitting region
            xi = abinitio[key].fit_range[1]
            xf = abinitio[key].fit_range[2]
            mask = (xi .< x_ai .< xf)
            #
            ## 
            x_ai = x_ai[mask]
            y_ai = y_ai[mask]
            #
            ## bra and ket labels
            # i = floor(Int64, key[2][1])
            # j = floor(Int64, key[2][2])
            #
            ## extract fitting flags, i.e. turn of parameter variation in fit
            p_excludeFromFit = fitting_object.fit
            #
            ## extract guesses for parameters
            p_guess = fitting_object.Rval
            p = deepcopy(p_guess)
            #
            ## determine the functional form to fit
            func = fitting_object.type
            #
            ## extract parameter bounds
            p_bounds = fitting_object.bounds
            #
            ## obtain parameters for the ComputePropert_viaParameters function
            Lval     = fitting_object.Lval
            units    = fitting_object.units
            sub_type = fitting_object.sub_type
            factor   = fitting_object.factor
            #
            ## if no parameters are to be fit then skip fitting step
            if (any(p_excludeFromFit .== 1))
                    #
                    ## perform optimization
                    o_ = optimize(p -> get_loss_general(x_ai, y_ai, Lval, units, key[1], sub_type, factor, func, p_excludeFromFit, p_guess, p_bounds, p), [p_guess...], options) 
                    p = Optim.minimizer(o_)
                    fitting_object.fitted_parameters = p
                    println(p)
            end
            #
            ## compute fitted curve
            x, fitted_object = ComputeProperty_viaParameters(r, func, Lval, p, key[1], units, sub_type, factor)
            #
            ## plot fitted curve
            plt.figure()
            plt.title("<"*string(i)*"|"*key[1]*"|"*string(j)*">")
            plt.plot(x_ai, y_ai, label = "abinitio")
            plt.plot(x, fitted_object, label = "fitted")
            # plt.xlim(x[1],x[end])
            plt.xlim(xi,xf)
            plt.legend()
        end
    end
end
#
function compute_spectroscopic_parameters(state, Rval; rmin = Calculation["grid"].range[1], rmax =Calculation["grid"].range[2], minima = false)
    #
    ## constants
    mu = 10.6613025642667
    amu = 1.660538921000E-24
    h = 6.626069570000E-27
    c = 29979245800
    #
    ## find true minimum
    if minima == false
        o_ = optimize(p -> ComputeProperty_viaParameters_SP(p..., 
                                                            Potential[state].type, 
                                                            Potential[state].Lval, 
                                                            Rval, 
                                                            Potential[state].obj_type, 
                                                            Potential[state].units, 
                                                            Potential[state].sub_type, 
                                                            Potential[state].factor,
                                                            state = state), 
                                                            rmin, 
                                                            rmax, 
                                                            [1.5],
                                                            Fminbox(LBFGS()))
        #
        Re = Optim.minimizer(o_)[1]
        Ve = Optim.minimum(o_)
    else
        Re, Ve = minima
    end
    #
    function centered_grid(x0, step, N)
        @assert N ≥ 1 "N must be at least 1"
        if N == 1
            return [x0]
        end
    
        # Make sure x0 is included exactly
        if isodd(N)
            # symmetric grid
            half = (N - 1) ÷ 2
            grid = x0 .+ step .* collect(-half:half)
        else
            # for even N, shift so x0 is one of the two middle points
            half = N ÷ 2
            grid = x0 .+ step .* collect(-(half-1):half)
        end
        return grid
    end
    #
    R_near_min = centered_grid(Re, 0.001, 16)
    #
    r_, V_ =  ComputeProperty_viaParameters(R_near_min, 
                                        Potential[state].type, 
                                        Potential[state].Lval, 
                                        Rval, 
                                        Potential[state].obj_type, 
                                        Potential[state].units, 
                                        Potential[state].sub_type, 
                                        Potential[state].factor)
    #
    ## compute derivatives up to 4th order
    der1 = FiniteDifferenceSP(R_near_min, V_, 8, 1)  #      -0.0000000001
    der2 = FiniteDifferenceSP(R_near_min, V_, 8, 2)  #  418798.9312634493
    der3 = FiniteDifferenceSP(R_near_min, V_, 8, 3)  # -2632204.1256825030
    der4 = FiniteDifferenceSP(R_near_min, V_, 8, 4)  #  13856329.20514625
    #
    ## compute Harmonic frequency
    we = sqrt(h * c * abs(der2) * 10^(16) / (mu * amu)) * (1/(2 * pi * c))
    #
    ## compute Rotational constant
    Be = h / (8 * pi * pi * mu * amu * Re * Re * 10^(-16) * c)
    #
    ## compute anharmonicity constant
    T1 = ((Be)^(2) * Re^(4)) / (12 * (we)^(5))
    #
    T2 = 10 * (Be) * (abs(der3))^(2) * Re^(2)
    #
    T3 = 3 * (abs(der4)) * we^(2)
    #
    xe = T1 * ( T2 - T3 )
    #
    r_, Ae =  ComputeProperty_viaParameters(1000, 
                                            Potential[state].type, 
                                            Potential[state].Lval, 
                                            Rval, 
                                            Potential[state].obj_type, 
                                            Potential[state].units, 
                                            Potential[state].sub_type, 
                                            Potential[state].factor)
    #
    return Ve, Re, Ae, Be, we, we*xe
end
#
function fit_spectroscopic_constants()
    function get_loss_general(state,
                              values, 
                              fixParams, 
                              p_init, 
                              p_bounds, 
                              p ; 
                              rmin=Calculation["grid"].range[1],
                              rmax= Calculation["grid"].range[2],
                              minima=false)
        #
        for (i,par) in enumerate(p)
                if fixParams[i] == 0.0
                        p[i] = p_init[i]
                else
                        p[i] = par 
                        #
                        if (p[i]<p_bounds[i][1])|(p[i]>p_bounds[i][2])
                                return 1e100
                        end
                end
        end
        #
        ## compute spectroscopic parameters
        Ve, Re, Ae, Be, we, wexe = compute_spectroscopic_parameters(state, p, rmin = rmin,rmax = rmax, minima=minima)
        #
        model = [Ve, Re, Ae, Be, we, wexe]
        #
        cost = 0
        #
        epsilon = 1e-6
        #
        for (idx,param) in enumerate(values)
            #
            if isnan(param)
                continue
            else
                cost += ((abs(param) - (model[idx]))/(abs(param) + abs(model[idx]) + epsilon))^2
            end
        end
        #
        # println("CAT",cost)
        return sqrt(cost)
    end
    #
    options = Optim.Options(show_trace = true)
    #
    for key in keys(SpectroscopicConstants)
        #
        ## extract fitting region
        rmin, rmax = SpectroscopicConstants[key].fit_range
        #
        if key in Calculation["method"].states
            fitting_object = Potential[key]
            #
            constants = SpectroscopicConstants[key].Lval
            values = SpectroscopicConstants[key].Rval
            #
            ## re order them based on output of constant calculator
            ordering = ["ve","re","ae","be","we","wexe"]
            #
            c = String[]
            v = Float64[]
            
            for constant in ordering
                i = findfirst(==(constant), constants)
                if isnothing(i)
                    push!(c, constant)
                    push!(v, NaN)
                else
                    push!(c, constant)
                    push!(v, values[i])
                end
            end
            #
            constants = c
            values    = v
            #
            ## extract fitting flags, i.e. turn of parameter variation in fit
            p_excludeFromFit = fitting_object.fit
            #
            ## extract guesses for parameters
            p_guess = fitting_object.Rval
            p = deepcopy(p_guess)
            #
            ## determine the functional form to fit
            func = fitting_object.type
            #
            ## extract parameter bounds
            p_bounds = fitting_object.bounds
            #
            ## obtain parameters for the ComputePropert_viaParameters function
            Lval     = fitting_object.Lval
            state  = fitting_object.ID
            #
            ## determine whether Re and Ve are parametrised
            min_param_funcs = ["morse_oscillator","emo"]
            if lowercase(fitting_object.type) in min_param_funcs
                minima = zeros(Float64,2)
                #
                for (idx, param) in enumerate(Lval)
                    if lowercase(param) == "re"
                        minima[1] = fitting_object.Rval[idx]
                    elseif lowercase(param) == "ve"
                            minima[2] = fitting_object.Rval[idx]
                    end
                end
            else
                minima = false
            end
            #
            ## if no parameters are to be fit then skip fitting step
            if (any(p_excludeFromFit .== 1))
                    #
                    ## perform optimization
                    o_ = optimize(p -> get_loss_general(state, values, p_excludeFromFit, p_guess, p_bounds, p, rmin=rmin, rmax=rmax, minima=minima), [p_guess...], options) 
                    p = Optim.minimizer(o_)
                    fitting_object.fitted_parameters = p
                    println(p)
            end
            #
            ## print fitted constants
            Ve, Re, Ae, Be, we, wexe = compute_spectroscopic_parameters(fitting_object.ID, p,rmin = rmin,rmax = rmax, minima=minima)
            #
            println()
            println("___ FITTED SPECTROSCOPIC PARAMETERS FOR STATE $state ___")
            println("Ve   = $Ve   cm-1      vs. ", v[1])
            println("Re   = $Re   Angstroms vs. ", v[2])
            println("Ae   = $Ae   cm-1      vs. ", v[3])
            println("Be   = $Be   cm-1      vs. ", v[4])
            println("we   = $we   cm-1      vs. ", v[5])
            println("wexe = $wexe cm-1      vs. ", v[6])
        end
    end
end




