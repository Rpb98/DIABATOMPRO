#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LIBRARIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
include("dependencies.jl")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ATOMIC MASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
include("atomic_masses.jl")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DEFINE OBJECT CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
include("Objects.jl")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
include("functions.jl")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INPUT READER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
include("input_reader.jl")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DIABATISER ROUTINES ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
include("Diabatiser.jl")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~ VIBRONIC SOLVER ROUTINES ~~~~~~~~~~~~~~~~~~~~~~~~~#
include("vibronic_solver.jl")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DIABATOMPRO FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
using PyPlot
#
## terminal print statement
print("▒█▀▀▄ ▀█▀ ░█▀▀█ ▒█▀▀█ █▀▀█ ▀▀█▀▀ █▀▀█ █▀▄▀█  ░░  ▒█▀▀█ █▀▀█ █▀▀█\n")    
print("▒█░▒█ ▒█░ ▒█▄▄█ ▒█▀▀▄ █▄▄█ ░░█░░ █░░█ █░▀░█  ▀▀  ▒█▄▄█ █▄▄▀ █░░█\n")    
print("▒█▄▄▀ ▄█▄ ▒█░▒█ ▒█▄▄█ ▀░░▀ ░░▀░░ ▀▀▀▀ ▀░░░▀  ░░  ▒█░░░ ▀░▀▀ ▀▀▀▀\n")    
print("\n")
print("𝚆𝚛𝚒𝚝𝚝𝚎𝚗 𝚋𝚢 𝚁𝚢𝚊𝚗 𝙱𝚛𝚊𝚍𝚢 (UCL) | 𝚟𝟷.𝟶 | 𝙻𝚊𝚜𝚝 𝚞𝚙𝚍𝚊𝚝𝚎𝚍 𝟶𝟸.𝟷𝟶.𝟸𝟶𝟸𝟹\n")
print("𝙴𝚖𝚊𝚒𝚕: 𝚛𝚢𝚊𝚗.𝚋𝚛𝚊𝚍y.𝟷𝟽@𝚞𝚌𝚕.𝚊𝚌.𝚞𝚔\n")
print("\n")
println("Conversion factors used by DIABATOM-PRO:")
println("              bohr  =       0.529177210920     Å")
println("                Eh  =   219474.63137080000  cm-1")
println("               erg  =   5.034117008194E+15  cm-1")
println("                eV  =      8065.5442959967  cm-1")
println("          kCal/mol  =       349.7550878997  cm-1")
println("            kJ/mol  =        83.5934722514  cm-1")
println("               THz  =        33.3564095198  cm-1")
println("a.u. of dipole ea0  =       2.541746363812 Debye")
print("\n")
#      
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUN INPUT READER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# fname =  "/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/DIABATOM_PRO_PACKAGE/github_dev/DIABATOMPRO/Supplementary/KH/KH.inp" #"/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/DIABATOM_PRO_PACKAGE/github_dev/DIABATOMPRO/Supplementary/CH/2Pi/CH_doublet_pi.inp" "/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/DIABATOM_PRO_PACKAGE/github_dev/DIABATOMPRO/Supplementary/SO/SO.inp"
read_file("/Users/ryanbrady/Documents/PostDoc/DIABATISATION/DIABATOM_PRO_PACKAGE/github_dev/DIABATOMPRO/src/tests/SO_d1Pi.inp") # <example_input_name> ../Supplementary/SO/SO_C_Cp_DMZ.inp ./TEST_input.inp"/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/DIABATOM_PRO_PACKAGE/github_dev/DIABATOMPRO/Supplementary/SO/SO_B.inp") ##../Supplementary/SO/SO_B.inp
#~~~~~~~~~~~~~~~~~~~~~~~~~ RUN HAMILTONIAN BUILDER ~~~~~~~~~~~~~~~~~~~~~~~~#
# include(joinpath(@__DIR__, "Build_Hamiltonian_Matrix.jl"))
include("Build_Hamiltonian_Matrix.jl")
#~~~~~~~~~~~~~~~~~~~~~~~ RUN DIABATISATION PIPELINES ~~~~~~~~~~~~~~~~~~~~~~#
#
## define the bond length vector
global r
r = LinRange(Calculation["grid"].range[1],
             Calculation["grid"].range[2],
             Calculation["grid"].npoints)  

r  = collect(r)

Ve, Re, Ae, Be, we, wexe = compute_spectroscopic_parameters(1, Potential[1].Rval, rmin = 1.35,rmax = 1.9, minima=false)


plt.figure()
plt.plot(r,PotMat[:,1,1])

v = [43902, 1.5303, 4.37878000000000E+04, 0.6752, 607.39, 45.37]
println()
println("___ FITTED SPECTROSCOPIC PARAMETERS FOR STATE 1 ___")
println("Ve   = $Ve   cm-1      vs. ", v[1])
println("Re   = $Re   Angstroms vs. ", v[2])
println("Ae   = $Ae   cm-1      vs. ", v[3])
println("Be   = $Be   cm-1      vs. ", v[4])
println("we   = $we   cm-1      vs. ", v[5])
println("wexe = $wexe cm-1      vs. ", v[6])

# compute vibronic energies and wavefunctions for a non-rotating molecule if
# the 'vibronic_solver' key is in Calculation
# run the diabatiser
# if Calculation["method"].abinitio_fit == true
#     fit_abinitio()
# elseif Calculation["method"].fit_spectroscopic_constants == true
#         fit_spectroscopic_constants()
# elseif (Calculation["method"].abinitio_fit == false) & (Calculation["method"].fit_spectroscopic_constants == false)
#     U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
   
#     fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])
#     plt.subplots_adjust(wspace=0, hspace=0)
#     for i=1:dim
#         if lowercase(Calculation["method"].regularisation) == "potential"
#             if i in Calculation["method"].states
#                 axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
#                 axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
#             end
#         else
#             for j=i:dim
#                 if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#                     axs[1,1].plot(r,Diabatic_Objects[Calculation["method"].regularisation][:,i,j],label="F"*string(i)*string(j))
#                     axs[1,1].plot(r,Objects[Calculation["method"].regularisation][:,i,j],"--")
#                     #
#                 end
#             end
#         end
#     end
#     #
#     for i=1:dim
#         for j=i+1:dim
#             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#                 axs[2,1].plot(r,Objects["nac"][:,i,j])
#             end
#         end
#     end


#     # plt.figure()
#     # plt.plot(r,map(x -> Objects["K_matrix"][x][1,1], collect(1:lastindex(r))))
#     # plt.xlabel("Bond Length, Angstroms")
#     # plt.ylabel("DBOC, cm-1")

#     # plt.figure()
#     # plt.plot(r,Diabatic_Objects["potential"][:,1,2])
#     # plt.xlabel("Bond Length, Angstroms")
#     # plt.ylabel("Diabatic (potential) coupling, cm-1")


#     # for i=1:dim
#     #     for j=i+1:dim
#     #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
                
#     #             # axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
#     #             axs[2,1].plot(r,Objects["nac"][:,i,j]) #,"--",alpha=0.5)
#     #         end
#     #     end
#     # end

#     # axs[2,1].plot(r,Objects["nac"][:,2,3],"k") #,"--",alpha=0.5)

#     # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
#     # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$")
#     # # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
#     # # axs[1,1].legend()
#     # # axs[2,1].legend()

#     # sf = Calculation["method"].states[end]
#     # Emax = Objects["potential"][end,sf,sf]

#     # si = Calculation["method"].states[1]
#     # Emin = minimum(Objects["potential"][:,si,si])

#     # # axs[1,1].set_ylim(Emin,1.1*Emax)
#     # axs[1,1].set_xlim(1.35,3.5)
# end




# save_diabatisation(Objects, Diabatic_Objects, lowercase(Calculation["method"].diabatisation), input_properties, "SO_19SaNadiab", special_name = "Cp")








# # # Ve, Re, Be, we, xe = compute_spectroscopic_parameters(1,rmin=1.2,rmax=2.0)

# # # println("Ve = $Ve cm-1, re = $Re Ang, we = $we cm-1, Be = $Be cm-1, xe = $xe")



# c = [Potential[1].Rval[15:end-4]...]
# plt.plot(r,repulsive.(r,c...))
# #
# plt.plot(r,EMO.(r,Potential[1].Rval[1:14]...))
# #
# plt.plot([1.5303],[43902],"x")
# #
# plt.ylim(40000,60000)
# plt.xlim(1.2,4)


# function d1Pi_params(p)
#     Ve,Re,Ae,B0,B1,A6,A7= p #,B2,B3,B4
#     #
#     Potential[1].Rval[1]  = Ve
#     Potential[1].Rval[2]  = Re
#     Potential[1].Rval[3]  = Ae
#     Potential[1].Rval[9]  = B0
#     Potential[1].Rval[10] = B1
#     Potential[1].Rval[22] = A6
#     Potential[1].Rval[23] = A7
#     # Potential[1].Rval[29] = gamma

#     #
#     x, y = ComputeProperty(Potential[1])
#     #
#     ## compute spectroscopic parameters
#     lower_bound = 1.2
#     upper_bound = 2.1
#     mask = (x .>= lower_bound) .& (x .<= upper_bound)
#     masked_E = y[mask]
#     #
#     V0_, min_index = findmin(masked_E)
#     #
#     ## now compute spectroscopic constants
#     Ve, Re, Be, we, xe = compute_spectroscopic_parameters(1,rmin = lower_bound,rmax = upper_bound)
#     #
#     ## compute barrier height
#     lower_bound = x[min_index]
#     upper_bound = x[end]
#     mask = (x .>= lower_bound) .& (x .<= upper_bound)
#     masked_E = y[mask]
#     Emax, max_index = findmax(masked_E)
#     #
#     barrier = Emax-Ve
#     #
#     ##
#     De = y[end]
#     #
#     model = [Ve,Re,barrier,we,xe,De]
#     #
#     return model
# end


# function cost_d1Pi(p,target)
#     model = d1Pi_params(p)
#     #
#     cost = 0
#     #
#     weights = [100,100,100,2,2,500]
#     #
#     for (idx,param) in enumerate(target)
#         cost += ((param - model[idx])/param)^2
#     end
#     #
#     return sqrt(cost)
# end

# function params(p)
#     #
#     Potential[2].Rval[3]  = p[1]
#     #
#     x, y = ComputeProperty(Potential[2])
#     #
#     ## compute spectroscopic parameters
#     lower_bound = 1.0
#     upper_bound = 4.0
#     mask = (x .>= lower_bound) .& (x .<= upper_bound)
#     masked_E = y[mask]
#     #
#     V0_, min_index = findmin(masked_E)
#     #
#     ## now compute spectroscopic constants
#     minparams = [Potential[2].Rval[2],Potential[2].Rval[1]]
#     Ve, Re, Be, we, xe = compute_spectroscopic_parameters(2, rmin = lower_bound,rmax = upper_bound, minima=minparams)
#     # println("CATS", Ve, Re, Be, we, xe)
#     #
#     ## compute barrier height
#     lower_bound = x[min_index]
#     upper_bound = x[end]
#     mask = (x .>= lower_bound) .& (x .<= upper_bound)
#     masked_E = y[mask]
#     Emax, max_index = findmax(masked_E)
#     #
#     barrier = Emax-Ve
#     #
#     ##
#     De = Potential[2].Rval[end]
#     #
#     model = [Ve,Re,we,we*xe,De]
#     #
#     return model
# end

# function cost_5Pi(p,target)
#     Ve,Re,we,wexe,De = params(p)
#     #
#     model = [we,wexe]
#     #
#     cost = 0
#     #
#     weights = [100,100,100,2,2,500]
#     #
#     for (idx,param) in enumerate(target)
#         cost += ((param - model[idx])/param)^2
#     end
#     #
#     return sqrt(cost)
# end

# 44963.0857
# 45540.45511

# p_guess = [Potential[2].Rval[3]]
   
# options = Optim.Options(show_trace = true)
# o_ = optimize(p -> cost_5Pi(p,[330.3, 12.26]), [p_guess...],options) 
# optimisedParameters = Optim.minimizer(o_)


# Potential[2].Rval[3]   = optimisedParameters[1]

# # Potential[1].Rval[13]   = optimisedParameters[8]

# expt = [43902, 1.5303, 1570, 607.39, 0.0747,4.37878000000000E+04]

# println(params(Potential[2].Rval[3]))

# # println(d1Pi_params(optimisedParameters).-expt)

# # #
# Potential[2].Rval = Potential[2].fitted_parameters
# x, y = ComputeProperty(Potential[2])

# plt.plot(x,y)
# plt.plot(r,PotMat[:,2,2])
# plt.ylim(40000,55000)
# plt.plot(abinitio[("poten",1)].Lval,abinitio[("poten",1)].Rval)







# p_guess = [Potential[1].Rval[1:3]..., Potential[1].Rval[9:10]...,Potential[1].Rval[22:23]...]
   
# options = Optim.Options(show_trace = true)
# o_ = optimize(p -> cost_d1Pi(p,[43902, 1.5303, 1570, 607.39, 0.0747,4.37878000000000E+04]), [p_guess...],options) 
# optimisedParameters = Optim.minimizer(o_)
# 
# optimisedParameters = [45599.1975116424,1.570262862330456,47092.63345541411,5.132817840445445,-5.272200632947819,-0.27016618510739876,-0.11121014641033768,-0.055224309749338354]

# optimisedParameters = [45528.61998581535,1.5683849423761331,46947.20328186522,5.219027550979376,-4.720131916979748,1.1121798995177096e6,0.45628387218938415]
# Potential[1].Rval[1]   = optimisedParameters[1]
# Potential[1].Rval[2]   = optimisedParameters[2]
# Potential[1].Rval[3]   = optimisedParameters[3]
# Potential[1].Rval[9]   = optimisedParameters[4]
# Potential[1].Rval[10]  = optimisedParameters[5]
# Potential[1].Rval[22]   = optimisedParameters[6]
# Potential[1].Rval[23]   = optimisedParameters[7]
# # Potential[1].Rval[13]   = optimisedParameters[8]

# expt = [43902, 1.5303, 1570, 607.39, 0.0747,4.37878000000000E+04]

# println(d1Pi_params(optimisedParameters))

# println(d1Pi_params(optimisedParameters).-expt)

# # #
# x, y = ComputeProperty(Potential[1])

# plt.plot(x,y)
# plt.plot(abinitio[("poten",1)].Lval,abinitio[("poten",1)].Rval)

# # data_matrix = hcat(x, y)

# # # 3. Write the matrix to the file using the specified delimiter
# # using DelimitedFiles
# # writedlm("./d1Pi_fitted2constants.dat", data_matrix, ' ')


# # # gamma = 0.6010394259406852 
# # # RC = 2.371569590950108
# # # plt.figure()
# # # plt.plot(r,lorentzian.(r,gamma,RC,1) .- duo_lorentzian.(r,gamma*2,RC,pi/2))
# # # plt.plot(r, duo_lorentzian.(r,gamma*2,RC,pi/2),"x")



# # # PotMat[:,1,1] = y

# # # lower_bound = 1.2
# # # upper_bound = 2.1
# # # mask = (r .>= lower_bound) .& (r .<= upper_bound)
# # # masked_E = PotMat[:,1,1][mask]
# # # Emin, min_index = findmin(masked_E)

# # # # println("COMPARISON")
# # # # println("Ve(expt.) = ",43902," Ve(model) = ",Emin, " with difference = ",43902 - Emin," UPDATED VE = ",Potential[1].Rval[1]+ 43902 - Emin)
# # # # println("Re(expt.) = ",1.5303," Re(model) = ",r[min_index], " with difference = ",1.5303 - r[min_index]," UPDATED RE = ",Potential[1].Rval[2]+ 1.5303 - r[min_index])

# # # # lower_bound = r[min_index]
# # # # upper_bound = r[end]
# # # # mask = (r .>= lower_bound) .& (r .<= upper_bound)
# # # # masked_E = PotMat[:,1,1][mask]
# # # # Emax, max_index = findmax(masked_E)

# # # # println("Barrier(expt.) >= ",1570," Barrier(model) = ",Emax - Emin, " with difference = ",1570 - (Emax - Emin))

# # # # println()

# # # plt.axhline(Emin+1570,color="k")



# # # data_matrix = hcat(r, PotMat[:,4,4])

# # # using DelimitedFiles

# # # writedlm("morphed_d1Pi.dat", data_matrix, ' ')




# # # compute vibronic wavefunctions
# # if haskey(Calculation, "vibronic_solver")
# #     #
# #     @time contr_vib_wfn, E_vib_contr = vibronic_eigensolver(collect(r), PotMat, Calculation["grid"].npoints, collect(keys(Potential)), Calculation["method"].atoms...)
# #     #
# #     plot_vibronic_solution(r, PotMat, [1], E_vib_contr,  contr_vib_wfn, [30])
# # end

# # # plt.figure()
# # # plt.plot(r,Diabatic_Objects["spin-orbit"][:,1,3])
# # # plt.plot(r,Diabatic_Objects["spin-orbit"][:,2,3])

# # # plt.plot(r,Objects["spin-orbit"][:,1,3],".")
# # # plt.plot(r,Objects["spin-orbit"][:,2,3],".")

# # # save_diabatisation(Objects, Diabatic_Objects, lowercase(Calculation["method"].diabatisation), input_properties, "SO_ai", special_name = "LSX_bC")


# # # compute v-dependent rotational constants and perturbation

# # function skew_lorentz(r, g0, r0, amp, a)
# #     #
# #     ## compute gamma
# #     if a != 0.0
# #         g = ref_sigmoid(r,g0,a,r0)
# #     else
# #         g = g0
# #     end
# #     #
# #     ## compute mixture
# #     return lorentzian(r, g, r0, amp)
# # end

# # function compute_alpha(r,U,f)
# #     alpha0 = zeros(length(r))
# #     #
# #     alpha_adi = [U[i,:,:] * [alpha0[i] 0 ; 0 f[i]] * U[i,:,:]' for i=1:lastindex(r)]
# #     alpha = [alpha_adi[i][1,1] for i=1:lastindex(r)]
# #     #
# #     return alpha
# # end

# # function B_cost(r, U, Tfac, contr_vib_wfn, vmax, Bv_exp, p)
# #     g0, r0, amp, a = p
# #     #
# #     ## compute functional correction
# #     f = skew_lorentz.(r, g0, r0, amp, a)
# #     #
# #     ## compute alpha
# #     alpha = compute_alpha(r,U,f)
# #     #
# #     ## compute rotational constants
# #     Bv     = zeros(vmax+1)
# #     alphav = zeros(vmax+1)
# #     #
# #     @views for v_idx=1:vmax+1
# #         #
# #         ## now for the rotational constant
# #         Bv[v_idx] =     Tfac * simps(contr_vib_wfn[i,v_idx,:] .* r.^(-2) .* contr_vib_wfn[i,v_idx,:],  r[1], r[end])
# #         alphav[v_idx] = Tfac * simps(contr_vib_wfn[i,v_idx,:] .* alpha .* r.^(-2) .* contr_vib_wfn[i,v_idx,:],  r[1], r[end])
# #         #
# #     end
# #     #
# #     ## compute residual
# #     res = Bv .+ alphav .- Bv_exp
# #     #
# #     ## compute RMSE
# #     return sqrt(sum(res.^2)/vmax)
# # end

# # function B_calc(r, U, Tfac, contr_vib_wfn, vmax, Bv_exp, p)
# #     g0, r0, amp, a = p
# #     #
# #     ## compute functional correction
# #     f = skew_lorentz.(r, g0, r0, amp, a)
# #     #
# #     ## compute alpha
# #     alpha = compute_alpha(r,U,f)
# #     #
# #     ## compute rotational constants
# #     Bv     = zeros(vmax+1)
# #     alphav = zeros(vmax+1)
# #     #
# #     @views for v_idx=1:vmax+1
# #         #
# #         ## now for the rotational constant
# #         Bv[v_idx] =     Tfac * simps(contr_vib_wfn[1,v_idx,:] .* r.^(-2) .* contr_vib_wfn[1,v_idx,:],  r[1], r[end])
# #         alphav[v_idx] = Tfac * simps(contr_vib_wfn[1,v_idx,:] .* alpha .* r.^(-2) .* contr_vib_wfn[1,v_idx,:],  r[1], r[end])
# #         #
# #     end
# #     #
# #     ## compute residual
# #     res = Bv .+ alphav .- Bv_exp
# #     #
# #     ## compute RMSE
# #     return sqrt(sum(res.^2)/vmax), res, Bv, alphav, alpha
# # end

# # # Tfac = KE_factor(Calculation["method"].atoms...)
# # # Bv_exp = [0.499129924662255,0.49399888067717795,0.48885188946379815,0.48371367391638587,0.4799,0.47379,0.47057,0.46442,0.459011,0.45298,0.44718,0.44233,0.4361,0.42904,0.4208,0.4102,0.38424,0.3281,0.30966,0.30345,0.29289,0.28425,0.27463,0.26306,0.25044,0.23779,0.22295,0.20496,0.18392,0.1621,0.1332]

# # # p_guess = [0.1, 3.0, -0.1, -10.0]


# # # options = Optim.Options(show_trace = true)
# # # o_ = optimize(p -> B_cost(r, U, Tfac, contr_vib_wfn, 30, Bv_exp,p), [p_guess...],options) 
# # # optimisedParameters = Optim.minimizer(o_)

# # # optimisedParameters = [0.19525109250083053,3.8457706650185166,-0.15779793910502107,3.818237863650165]

# # # cost, res, Bv, alphav, alpha = B_calc(r, U, Tfac, contr_vib_wfn, 30, Bv_exp,optimisedParameters)


# # # open("B_bob-rot.dat", "w") do io
# # #     # Loop through the vectors and write each pair of values
# # #     for i in 1:length(r)
# # #         # Write the x and y values, separated by a space, to a new line
# # #         println(io, r[i], " ", alpha[i])
# # #     end
# # # end

# # # # plt.figure()
# # # # plt.plot(r,alpha)

# # # println("COST = ", cost)


# # # vmax = Calculation["vibronic_solver"].contraction
# # # contracted_vibronic_dim = sum(vmax)

# # # i = 1

# # # Bv     = zeros(31)
# # # alphav = zeros(31)

# # # Tfac = KE_factor(Calculation["method"].atoms...)
# # # @views for v_idx=1:31
# # #     #
# # #     ## now for the rotational constant
# # #     Bv[v_idx] =     Tfac * simps(contr_vib_wfn[i,v_idx,:] .* r.^(-2) .* contr_vib_wfn[i,v_idx,:],  r[1], r[end])
# # #     alphav[v_idx] = Tfac * simps(contr_vib_wfn[i,v_idx,:] .* alpha .* r.^(-2) .* contr_vib_wfn[i,v_idx,:],  r[1], r[end])
# # #     #
# # # end


# # # J = 50
# # # v = collect(1:31)
# # # plt.figure()
# # # # plt.plot(v,Bv .+ alphav,"bx")
# # # # plt.plot(v,Bv_exp ,"r.")

# # # plt.plot(v,(Bv .+ alphav),"bx", label="Theory")
# # # # plt.plot(v,(Bv .- Bv_exp)           .* J*(J+1),"r.")
# # # plt.plot(v,Bv_exp,"r.",label="Expt.")
# # # # plt.ylim(0,0.6)
# # # plt.xlim(3,32)
# # # plt.xlabel("v")
# # # plt.ylabel(L"$B_v$, cm$^{-1}$")
# # # plt.legend()
# # # plt.title("With BOB-ROT correction")

# # # plt.savefig("Bv_with_bob-rot.png",dpi=300,bbox_inches="tight")


# # # plt.figure()
# # # plt.plot(v,Bv,"bx", label="Theory")
# # # plt.plot(v,Bv_exp,"r.",label="Expt.")
# # # plt.xlim(3,32)
# # # plt.xlabel("v")
# # # plt.ylabel(L"$B_v$, cm$^{-1}$")
# # # plt.legend()
# # # plt.title("Without BOB-ROT correction")

# # # plt.savefig("Bv_without_bob-rot.png",dpi=300,bbox_inches="tight")

# # # plt.plot(v,alphav)


# # # M = zeros(Float64,31, length(r))

# # # for v=1:31
# # #     for geom=1:lastindex(r)
# # #         M[v,geom] = contr_vib_wfn[i,v,geom]^2
# # #     end
# # # end

# # # F = svd(M)
# # # U = F.U
# # # S = F.S
# # # Vt = F.Vt
# # # V = collect(Vt') # Get the non-transposed V matrix


# # # # Step 3: Compute the pseudoinverse of M
# # # # Julia has a built-in function 'pinv' which is the most direct way.
# # # # It handles the truncation of small singular values automatically.
# # # M_pseudoinverse = pinv(M)


# # # # Step 4: Solve for the unknown vector 'a' using the pseudoinverse
# # # # The solution a is the minimum-norm solution that fits the data.
# # # # This 'a' vector contains the values of alpha(r) at each grid point.
# # # b = Bv .- Bv_exp
# # # a = M_pseudoinverse * b

# # # plt.figure()
# # # plt.plot(r,a)






# # # #
# # # ## compute expectation values: Psi_v = psi[state,v_idx,:]
# # # expec = []
# # # for v=1:31
# # #     Psi_v = contr_vib_wfn[1,v,:]
# # #     #
# # #     integrand = Psi_v .* r .* Psi_v
# # #     #
# # #     ## integrate
# # #     rexpec = simps(integrand, r[1], r[end])
# # #     # println(rexpec)
# # #     #
# # #     println("<$(v-1)|r|$(v-1)> =", rexpec)
# # # end

# # # diffs = [0.1,1,2.5,5,10,20,50,100]

# # # regions = LinRange(0.01,80000,80000)

# # # plt.figure()

# # # for diff in diffs
# # #     power = regions/diff
# # #     plt.plot(regions,power,label=diff)
# # # end
# # # plt.yscale("log")

# # # plt.legend()


# # ########################## ROVIBRONIC CODE #####################################


# # # using Profile

# # # state = 1
# # # J = 2.0
# # # #
# # # ## generate absolute AM Hunds case (a) basis & reflection operator phases
# # # Lambda, Sigma, Omega, S, phase = generate_allowed_AM_values_2(state, J)
# # # #
# # # for Tau in ["+","-"]
# # #     for i=1:lastindex(Lambda)
# # #         #
# # #         ## compute symmetrised basis
# # #         println("1/√2 (|$(Lambda[i]),$(Sigma[i]),$(Omega[i])> $(Tau) $(phase[i])*|$(-Lambda[i]),$(-Sigma[i]),$(-Omega[i])>)")
# # #     end
# # # end

# # # Omega = -2,-1,0,0,1,2
# # # Jz = [-2 0 0 0 0 0 ; 0 -1 0 0 0 0 ; 0 0 0 0 0 0; 0 0 0 0 0 0 ; 0 0 0 0 1 0 ; 0 0 0 0 0 2]

# # # U = 1/sqrt(2) * [1 0 0 0 0 -1 ; 0 1 0 0 -1 0; 0 0 sqrt(2) 0 0 0 ; 0 0 0 sqrt(2) 0 0 ; 0 1 0 0 1 0 ; 1 0 0 0 0 1]

# # # Jz_trans = U * Jz * U'

# # function generate_symmetrized_AM_basis(state::Int64, J::Float64)
# #     #
# #     ## Get raw Hunds case (a) basis
# #     Lambda, Sigma, Omega, S = generate_allowed_AM_values(state, J)
# #     #
# #     ## compute inversion QN
# #     if (Lambda == 0)&(Potential[state].symmetry == "-")
# #         s=1
# #     else
# #         s=0
# #     end
# #     #
# #     N = length(Lambda)
# #     used = falses(N)
# #     #
# #     sym_basis = []
# #     #
# #     for i in 1:N
# #         if used[i]
# #             continue
# #         end
# #         # Look for parity partner j
# #         for j in (i+1):N
# #             if !used[j] &&
# #                 isapprox(Lambda[i] + Lambda[j], 0.0) &&
# #                 isapprox(Sigma[i] + Sigma[j], 0.0) &&
# #                 isapprox(Omega[i] + Omega[j], 0.0)
# #                 #
# #                 ## Compute ε phase
# #                 exponent = Int(round(s - Lambda[i] + S - Sigma[i] + J - Omega[i]))
# #                 ε = (-1)^exponent
# #                 #
# #                 #
# #                 ## Store ± combinations as symbolic expressions (or any structure you prefer)
# #                 push!(sym_basis, ("+", (Lambda[i], Sigma[i], Omega[i]), (Lambda[j], Sigma[j], Omega[j]), ε))
# #                 push!(sym_basis, ("-", (Lambda[i], Sigma[i], Omega[i]), (Lambda[j], Sigma[j], Omega[j]), ε))
# #                 #
# #                 used[i] = true
# #                 used[j] = true
# #                 break
# #             end
# #         end
# #     end
# #     #
# #     return sym_basis
# # end

# # function compute_symmetrised_rovibronic_energies(J,Tau)
# #     QN_book, QN_index_map = quantum_number_bookkeeping(J,Tau)

# #     contr_vib_wfn, E_vib_contr = vibronic_eigensolver(collect(r), PotMat, Calculation["grid"].npoints, collect(keys(Potential)), Calculation["method"].atoms...)

# #     T, V, B = build_coupled_vibronic_Hamiltonian(contr_vib_wfn)

# #     Hrot_list = build_rotational_parity_subBlocks(J,Tau)

# #     Hrot, spin_rot_dims, tot_rovibronic_dim = build_Hrot(B, Hrot_list)

# #     Hvibronic =  build_Hvibronic_embedding(T, V, spin_rot_dims, tot_rovibronic_dim)

# #     Htot = Hvibronic + Hrot

# #     E, eigvec = eigen(Htot)

# #     #

# #     return QN_book, QN_index_map, E, eigvec, contr_vib_wfn, E_vib_contr
# # end


# # # Jlist = [0.0]

# # # QN_books = []
# # # QN_index_maps = []
# # # eners = []

# # # minE = []

# # # points = LinRange(100,1010,5)

# # # v1 = []

# # # for J in Jlist
# # #     for Tau in [1]

# # #     # QN_book, QN_index_map = quantum_number_bookkeeping(J,Tau)

# # #     # contr_vib_wfn, E_vib_contr = vibronic_eigensolver(collect(r), PotMat, Calculation["grid"].npoints, collect(keys(Potential)), Calculation["method"].atoms...)

# # #     # T, V, B = build_coupled_vibronic_Hamiltonian(contr_vib_wfn)

# # #     # Hrot_list = build_rotational_parity_subBlocks(J,Tau)

# # #     # Hrot, spin_rot_dims, tot_rovibronic_dim = build_Hrot(B, Hrot_list)

# # #     # Hvibronic =  build_Hvibronic_embedding(T, V, spin_rot_dims, tot_rovibronic_dim)

# # #     # Htot = Hvibronic + Hrot

# # #     # E, eigvec = eigen(Htot)

# # #     # if J == 0.0
# # #     #     ZPE = minimum(E)
# # #     # end
# # #         # Calculation["grid"].npoints = N
# # #         # r = LinRange(Calculation["grid"].range[1],
# # #         #              Calculation["grid"].range[2],
# # #         #              Calculation["grid"].npoints)  
# # #         # include("Build_Hamiltonian_Matrix.jl")
# # #         QN_book, QN_index_map, E, eigvec, contr_vib_wfn, E_vib_contr = compute_symmetrised_rovibronic_energies(J,Tau)

# # #         push!(QN_books, QN_book)
# # #         push!(QN_index_maps, QN_index_map)
# # #         push!(eners, E)
# # #         push!(minE, minimum(E))
# # #         push!(v1,E[2]-minimum(E))
# # #         print("Ev=1 = ",E[2]-minimum(E))
    
# # #     end
# # # end

# # # print(eners)
# # # ZPE = minimum(minE)

# # # eners = [E .- ZPE for E in eners]



# # #
# # # open("J0_SO.eners", "w") do io
# # #     # println(io, "|--------------------------------------------------------------------------|")
# # #     println(io, "    n   J        E         state    v    S    Lambda    Sigma    Omega    Parity ")
# # #     # println(io, "|--------------------------------------------------------------------------|")
# # #     #
# # #     for (idx,J) in enumerate(Jlist)
# # #         E = eners[idx]
# # #         QN_book = QN_books[idx]
# # #         QN_index_map = QN_index_maps[idx]
# # #         for n in 1:lastindex(E)
# # #             qn = QN_book[n]
# # #             state   = qn[QN_index_map["state"]]
# # #             v       = qn[QN_index_map["v"]]
# # #             Lambda  = qn[QN_index_map["Lambda"]]
# # #             Sigma   = qn[QN_index_map["Sigma"]]
# # #             J       = qn[QN_index_map["J"]]
# # #             Omega   = qn[QN_index_map["Omega"]]
# # #             Parity  = qn[QN_index_map["Parity"]]
# # #             S       = (Potential[state].mult - 1) / 2

# # #             if Parity == 1
# # #                 Tau = "+"
# # #             elseif Parity == -1
# # #                 Tau = "-"
# # #             end

# # #             # Tau = "X"

# # #             @printf(io, "  %3d  %2.1f  %10.6f   %5d  %3d   %2.1f   %6d   %6.1f   %6.1f   %3s  \n",
# # #                     n, J, E[n], state, v, S, Lambda, Sigma, Omega, Tau)
# # #         end
# # #     end
# # # end

# # #
# # # H = T + V

# # # eig = eigen(H)


# # # U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))

# # # plt.figure()
# # # for key in keys(Dipole)
# # #     i, j = key
# # #     #
# # #     plt.plot(r,Diabatic_Objects["dipole"][:,i,j])
# # #     #
# # #     plt.plot(r,Objects["dipole"][:,i,j],"--")
# # # end
# #     #


# # # plt.ylim(40000,60000)

# # # plt.figure()
# # # plt.plot(r, PotMat[:,3,3])
# # # plt.plot(abinitio[("poten",3)].Lval, abinitio[("poten",3)].Rval,"--")
# # # # plt.plot(r, PotMat[:,4,4])

# # # # plt.plot(r, PotMat[:,1,1],"--")
# # # # plt.plot(r, PotMat[:,2,2],"--")

# # # # plt.plot(r, PotMat[:,5,5])
# # # # plt.plot(r, PotMat[:,6,6])
# # # # plt.plot(r, PotMat[:,7,7])

# # # plt.ylim(35000,75000)






# # # #
# # # ## initialise diabatom object
# # # diabatom = Dict()
# # # diabatom["r"] = r
# # # diabatom["Hamiltonian"] = Hamiltonian
# # # diabatom["Objects"] = Objects
# # # diabatom["Calculation"] = Calculation
# # # diabatom["Potential"] = Potential
# # # diabatom["SpinOrbit"] = SpinOrbit
# # # diabatom["EAMC"] = EAMC
# # # diabatom["Dipole"] = Dipole
# # # diabatom["NonAdiabaticCoupling"] = NonAdiabaticCoupling

# # # run the diabatiser
# # # if Calculation["method"].abinitio_fit == true
# # #     fit_abinitio()
# # # elseif Calculation["method"].abinitio_fit == false
# # #     U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
    
# #     # fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])

# #     # plt.subplots_adjust(wspace=0, hspace=0)

# #     # # axs[1,1].set_title("MRCI aug-cc-pVQZ-X2C | occ = 8330, closed = 5110")
# #     # for i=1:dim
# #     #     if lowercase(Calculation["method"].regularisation) == "potential"
# #     #         if i in Calculation["method"].states
# #     #             axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
# #     #             axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
# #     #         end
# #     #     else
# #     #         for j=i:dim
# #     #             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# #     #                 axs[1,1].plot(r,Diabatic_Objects[Calculation["method"].regularisation][:,i,j],label="F"*string(i)*string(j))
# #     #                 axs[1,1].plot(r,Objects[Calculation["method"].regularisation][:,i,j],"--")
# #     #             end
# #     #         end
# #     #     end
# #     # end

# #     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], color="green",label= L"$V^{\rm (d)}_1$")
# #     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], color="orange",label=L"$V^{\rm (d)}_2$")
# #     # # axs[1,1].plot(r,Objects["potential"][:,2,2],color="red","--")
# #     # # axs[1,1].plot(r,Objects["potential"][:,3,3],color="blue","--")

# #     # # inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25])  # [left, bottom, width, height]

# #     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,2,2], color="green")
# #     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,3,3], color="orange")
# #     # # inset_ax.plot(r,Objects["potential"][:,2,2],color="red","--")
# #     # # inset_ax.plot(r,Objects["potential"][:,3,3],color="blue","--")

# #     # # inset_ax.set_xlim(1.8, 2.1)  # Zoom into x-range
# #     # # inset_ax.set_ylim(50000,63000)  # Zoom into y-range

# #     # # axs[1,1].plot([1.8,1.8],[5e4,6.3e4], color="grey"     ,alpha=0.5)
# #     # # axs[1,1].plot([2.1,2.1],[5e4,6.3e4], color="grey"     ,alpha=0.5)
# #     # # axs[1,1].plot([1.8,2.1],[5e4,5e4], color="grey"       ,alpha=0.5)
# #     # # axs[1,1].plot([1.8,2.1],[6.3e4,6.3e4], color="grey"   ,alpha=0.5)
# #     # # axs[1,1].plot([2.1,2.668],[5e4,5e4], color="grey"     ,alpha=0.5)
# #     # # axs[1,1].plot([2.1,2.668],[6.3e4,7.09e4], color="grey",alpha=0.5)

# #     # # inset_ax.set_xticklabels([])  # Remove x-axis numbers
# #     # # inset_ax.set_yticklabels([])  # Remove y-axis numbers

# #     # # # axs[1,1].set_xlabel("Bond Length")
# #     # # axs[1,1].set_ylabel(L"Potential Energy, $cm^{-1}$")

# #     # for i=1:dim
# #     #     for j=i+1:dim
# #     #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# #     #             # axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
# #     #             axs[2,1].plot(r,Objects["nac"][:,i,j]) #,"--",alpha=0.5)
# #     #         end
# #     #     end
# #     # end

# #     # # axs[2,1].plot(r,Objects["nac"][:,2,3],"k") #,"--",alpha=0.5)


# #     # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
# #     # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$")
# #     # # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
# #     # # axs[1,1].legend()
# #     # # axs[2,1].legend()

# #     # sf = Calculation["method"].states[end]
# #     # Emax = Objects["potential"][end,sf,sf]

# #     # si = Calculation["method"].states[1]
# #     # Emin = minimum(Objects["potential"][:,si,si])

# #     # axs[1,1].set_ylim(Emin,1.1*Emax)
# #     # axs[1,1].set_xlim(1.35,3.5)

# #     # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/property_based_pot_guess.png",dpi=300,bbox_inches="tight")
# # # end




# # # fig, axs = plt.subplots(2,2,sharex=true,figsize=[5,3])
# # # plt.subplots_adjust(wspace=0.4, hspace=0)
# # # for i=1:dim
# # #     axs[1,1].plot(r,Objects["potential"][:,i,i])
# # #     axs[1,2].plot(r,Objects["potential"][:,i,i],alpha=0.25)

# # #     #
# # #     axs[1,2].plot(r,Diabatic_Objects["potential"][:,i,i])
# # #     #
# # #     for j=i+1:dim
# # #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #             axs[2,1].plot(r,Objects["nac"][:,i,j], alpha=0.5)
# # #             # axs[2,1].plot(r,[Objects["regularised_nac"][idx][i,j] for idx=1:lastindex(r)])
# # #             axs[2,2].plot(r,Diabatic_Objects["potential"][:,i,j])
# # #         end
# # #     end
# # # end

# # # col = ["k","red","orange","salmon","brown","blue","green"]

# # # col2 = ["blue","green","orange"]
# # # fig, axs = plt.subplots(2,1,sharex=true,figsize=[5,4])
# # # plt.subplots_adjust(wspace=0.4, hspace=0)

# # # global count,count2
# # # count = 0
# # # count2 = 0
# # # for i=1:dim
# # #     global count, count2
# # #     count+=1
# # #     axs[1,1].plot(r,Objects["potential"][:,i,i],linewidth=1,color=col[count],label=Potential[i].name)
# # #     #
# # #     for j=i+1:dim
# # #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #             global count2
# # #             count2+=1
# # #             axs[2,1].plot(r,Objects["nac"][:,i,j],color=col2[count2],label=string(i)*string(j))

# # #             Wij_ai = abinitio[("NAC",[i,j])]

# # #             idx = closest_value_index(Wij_ai.Lval, 2.8)

# # #             axs[2,1].plot(Wij_ai.Lval[idx:end],Wij_ai.Rval[idx:end],"--",color="k",alpha=1)
# # #         end
# # #     end
# # # end
# # # # axs[1,1].legend()
# # # axs[1,1].set_ylim(0,45000)
# # # axs[1,1].set_xlim(r[1],16)
# # # axs[1,1].set_ylabel(L"Potential Energy, cm$^{-1}$")
# # # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")


# # # axs[1,1].text(3.85,6500,L"$X^1\Sigma^+$",fontsize=8)
# # # axs[1,1].text(2.5,11000,L"$(1)^3\Sigma^+$",color="salmon",fontsize=8)
# # # axs[1,1].text(5.2,18000,L"$(2)^1\Sigma^+$",color="red",fontsize=8)
# # # axs[1,1].text(4.2,22500,L"$(1)^3\Pi$",color="green",fontsize=8)
# # # axs[1,1].text(2.6,22000,L"$(1)^1\Pi$",color="blue",fontsize=8)
# # # axs[1,1].text(3.5,33500,L"$(2)^3\Sigma^+$",color="brown",fontsize=8)
# # # axs[1,1].text(6.5,33500,L"$(3)^1\Sigma^+$",color="orange",fontsize=8)

# # # axs[1,1].text(13,15500,"K(4s)+H")
# # # axs[1,1].text(13,27500,"K(4p)+H")
# # # axs[1,1].text(13,36500,L"K$^+$+H$^-$",rotation=3)

# # # # axs[2,1].legend()

# # # axs[2,1].text(6,-0.2,L"$\langle X^1\Sigma^+|\frac{d}{dr}|(2)^1\Sigma^+\rangle$",color="blue",fontsize=9)
# # # axs[2,1].text(9.6,0.2,L"$\langle (2)^1\Sigma^+|\frac{d}{dr}|(3)^1\Sigma^+\rangle$",color="green",fontsize=9)
# # # axs[2,1].text(1.4,0.1,L"$\langle X^1\Sigma^+|\frac{d}{dr}|(3)^1\Sigma^+\rangle$",color="orange",fontsize=9)

# # # # plt.savefig("./KH_abinitio_PECs_NACs.png",bbox_inches="tight",dpi=300)




# # # lab_dict = Dict{Int, String}(
# # #       1 => L"$X^1\Sigma^+$",
# # #       2 => L"$(2)^1\Sigma^+$",
# # #       3 => L"$(3)^1\Sigma^+$",
# # #       4 => L"$(1)^3\Sigma^+$",
# # #       5 => L"$(2)^3\Sigma^+$",
# # #       6 => L"$(1)^1\Pi$",
# # #       7 => L"$(1)^3\Pi$"
# # #   )
# # # plt.figure(figsize=[5,3])
# # # for key in keys(SpinOrbit)
# # #     i, j = key
# # #     #
# # #     if i==j
# # #         sym = L"$_{\rm z}$"
# # #     else
# # #         sym = L"$_{\rm x}$"
# # #     end
# # #     #
# # #     plt.plot(r,Objects["spin-orbit"][:,i,j],label=L"\langle"*lab_dict[i]*L"|SO"*sym*"|"*lab_dict[j]*L"$\rangle$")
# # # end
# # # plt.legend(loc="center right",fontsize=8)

# # # plt.xlabel(L"Bond Length, $\rm \AA$")
# # # plt.ylabel(L"Spin-Orbit Coupling, cm$^{-1}$")

# # # plt.savefig("../Supplementary/KH/KH_abinitio_SOCs.png",bbox_inches="tight",dpi=300)


# # # plt.figure()
# # # for key in keys(EAMC)
# # #     i, j = key
# # #     #
# # #     plt.plot(r,Objects["lx"][:,i,j])
# # # end


# # # plt.figure(figsize=[5,3])
# # # for key in keys(EAMC)
# # #     i, j = key
# # #     #
# # #     sym = L"$_{\rm x}$"
# # #     #
# # #     plt.plot(r,Objects["lx"][:,i,j],label=L"\langle"*lab_dict[i]*L"|L"*sym*"|"*lab_dict[j]*L"$\rangle$")
# # # end
# # # plt.legend(loc="center right",fontsize=8)

# # # plt.xlabel(L"Bond Length, $\rm \AA$")
# # # plt.ylabel(L"Electronic Angular Momentum, $\hbar$")

# # # plt.savefig("../Supplementary/KH/KH_abinitio_LX.png",bbox_inches="tight",dpi=300)






# # # fig, axs = plt.subplots(2,1,figsize=[5,6],gridspec_kw=Dict("height_ratios" => [1.5, 1]))
# # # plt.subplots_adjust(hspace=0)
# # # for key in keys(Dipole)
# # #     i, j = key
# # #     #
# # #     name = Dipole[key].name
# # #     #
# # #     if occursin("dmx",lowercase(name))
# # #         C = L"$\mu_{\rm x}$"
# # #         axs[2,1].plot(r,Objects["dipole"][:,i,j],label=L"\langle"*lab_dict[i]*L"|"*C*"|"*lab_dict[j]*L"$\rangle$")

# # #     elseif occursin("dmz",lowercase(name))
# # #         C = L"$\mu_{\rm z}$"
# # #         axs[1,1].plot(r,Objects["dipole"][:,i,j],label=L"\langle"*lab_dict[i]*L"|"*C*"|"*lab_dict[j]*L"$\rangle$")

# # #     end
# # #     #
# # # end
# # # axs[1,1].legend(loc="center right",fontsize=8)
# # # axs[2,1].legend(loc="center right",fontsize=8)

# # # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
# # # axs[2,1].set_ylabel("Dipole Moment, Debye")
# # # axs[1,1].set_ylabel("Dipole Moment, Debye")


# # # plt.savefig("../Supplementary/KH/KH_abinitio_DMCs.png",bbox_inches="tight",dpi=300)





# # # col2 = ["blue","green","orange"]
# # # plt.figure()
# # # global count,count2
# # # count2 = 0
# # # for i=1:dim
# # #     global count2
# # #     #
# # #     for j=i+1:dim
# # #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #             global count2
# # #             count2+=1
# # #             #
# # #             Wij_ai = abinitio[("NAC",[i,j])]
# # #             idx = closest_value_index(Wij_ai.Lval, 3.5)
# # #             idxf = closest_value_index(Wij_ai.Lval, 9)

# # #             #
# # #             spl = Spline1D(r,Objects["nac"][:,i,j])
# # #             W_func = spl(Wij_ai.Lval)
# # #             #
# # #             # plt.plot(r,Objects["nac"][:,i,j],color=col2[count2],label=string(i)*string(j))

# # #             diff = (Wij_ai.Rval[idx:idxf] .- W_func[idx:idxf])
# # #             #
# # #             RMS = sqrt(sum(diff.^2)/length(diff))
# # #             println("RMS for ",i," ",j," = ",RMS)
# # #             plt.plot(Wij_ai.Lval[idx:idxf],(Wij_ai.Rval[idx:idxf] .- W_func[idx:idxf]),color=col2[count2],alpha=1)
# # #         end
# # #     end
# # # end






# # # fig,axs = plt.subplots(2,2,sharex=true)
# # # axs[1,1].plot(r,Objects["potential"][:,2,2])
# # # axs[1,1].plot(r,Objects["potential"][:,3,3])
# # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2],color="k",alpha=0.25)
# # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3],color="k",alpha=0.25)
# # # # plt.ylim(0.9*minimum(Objects["potential"][:,2,2]),1.1*Objects["potential"][end,3,3])
# # # axs[1,1].set_ylim(50000,60000)
# # # axs[1,1].set_xlim(1.4,2.5)

# # # axs[2,1].plot(r,Diabatic_Objects["potential"][:,2,3],color="k",alpha=0.25)




# # # NonAdiabaticCoupling[[2,3]].Rval[2] = 1.90898597152284660

# # # if Calculation["method"].abinitio_fit == true
# # #     fit_abinitio()
# # # elseif Calculation["method"].abinitio_fit == false
# # #     U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
# # #     #
# # #     # fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])

# # #     # plt.subplots_adjust(wspace=0, hspace=0)

# # #     # # axs[1,1].set_title("MRCI aug-cc-pVQZ-X2C | occ = 8330, closed = 5110")
# # #     # for i=1:dim
# # #     #     if lowercase(Calculation["method"].regularisation) == "potential"
# # #     #         if i in Calculation["method"].states
# # #     #             axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
# # #     #             axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
# # #     #         end
# # #     #     else
# # #     #         for j=i:dim
# # #     #             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #     #                 axs[1,1].plot(r,Diabatic_Objects[Calculation["method"].regularisation][:,i,j],label="F"*string(i)*string(j))
# # #     #                 axs[1,1].plot(r,Objects[Calculation["method"].regularisation][:,i,j],"--")
# # #     #             end
# # #     #         end
# # #     #     end
# # #     # end

# # #     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], color="green",label= L"$V^{\rm (d)}_1$")
# # #     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], color="orange",label=L"$V^{\rm (d)}_2$")
# # #     # # axs[1,1].plot(r,Objects["potential"][:,2,2],color="red","--")
# # #     # # axs[1,1].plot(r,Objects["potential"][:,3,3],color="blue","--")

# # #     # # inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25])  # [left, bottom, width, height]

# # #     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,2,2], color="green")
# # #     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,3,3], color="orange")
# # #     # # inset_ax.plot(r,Objects["potential"][:,2,2],color="red","--")
# # #     # # inset_ax.plot(r,Objects["potential"][:,3,3],color="blue","--")

# # #     # # inset_ax.set_xlim(1.8, 2.1)  # Zoom into x-range
# # #     # # inset_ax.set_ylim(50000,63000)  # Zoom into y-range

# # #     # # axs[1,1].plot([1.8,1.8],[5e4,6.3e4], color="grey"     ,alpha=0.5)
# # #     # # axs[1,1].plot([2.1,2.1],[5e4,6.3e4], color="grey"     ,alpha=0.5)
# # #     # # axs[1,1].plot([1.8,2.1],[5e4,5e4], color="grey"       ,alpha=0.5)
# # #     # # axs[1,1].plot([1.8,2.1],[6.3e4,6.3e4], color="grey"   ,alpha=0.5)
# # #     # # axs[1,1].plot([2.1,2.668],[5e4,5e4], color="grey"     ,alpha=0.5)
# # #     # # axs[1,1].plot([2.1,2.668],[6.3e4,7.09e4], color="grey",alpha=0.5)

# # #     # # inset_ax.set_xticklabels([])  # Remove x-axis numbers
# # #     # # inset_ax.set_yticklabels([])  # Remove y-axis numbers

# # #     # # # axs[1,1].set_xlabel("Bond Length")
# # #     # # axs[1,1].set_ylabel(L"Potential Energy, $cm^{-1}$")

# # #     # for i=1:dim
# # #     #     for j=i+1:dim
# # #     #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #     #             # axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
# # #     #             axs[2,1].plot(r,Objects["nac"][:,i,j]) #,"--",alpha=0.5)
# # #     #         end
# # #     #     end
# # #     # end

# # #     # # axs[2,1].plot(r,Objects["nac"][:,2,3],"k") #,"--",alpha=0.5)


# # #     # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
# # #     # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$")
# # #     # # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
# # #     # # axs[1,1].legend()
# # #     # # axs[2,1].legend()

# # #     # sf = Calculation["method"].states[end]
# # #     # Emax = Objects["potential"][end,sf,sf]

# # #     # si = Calculation["method"].states[1]
# # #     # Emin = minimum(Objects["potential"][:,si,si])

# # #     # axs[1,1].set_ylim(Emin,1.1*Emax)
# # #     # # axs[1,1].set_xlim(1.35,3.5)

# # #     # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/property_based_pot_guess.png",dpi=300,bbox_inches="tight")
# # # end


# # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2])
# # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3])

# # # axs[2,1].plot(r,Diabatic_Objects["potential"][:,2,3])







# # # Calculation["method"].states = [1,2]


# # # if Calculation["method"].abinitio_fit == true
# # #     fit_abinitio()
# # # elseif Calculation["method"].abinitio_fit == false
# # #     U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
# # #     #
# # #     # fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])

# # #     # plt.subplots_adjust(wspace=0, hspace=0)

# # #     # # axs[1,1].set_title("MRCI aug-cc-pVQZ-X2C | occ = 8330, closed = 5110")
# # #     # for i=1:dim
# # #     #     if lowercase(Calculation["method"].regularisation) == "potential"
# # #     #         if i in Calculation["method"].states
# # #     #             axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
# # #     #             axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
# # #     #         end
# # #     #     else
# # #     #         for j=i:dim
# # #     #             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #     #                 axs[1,1].plot(r,Diabatic_Objects[Calculation["method"].regularisation][:,i,j],label="F"*string(i)*string(j))
# # #     #                 axs[1,1].plot(r,Objects[Calculation["method"].regularisation][:,i,j],"--")
# # #     #             end
# # #     #         end
# # #     #     end
# # #     # end

# # #     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], color="green",label= L"$V^{\rm (d)}_1$")
# # #     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], color="orange",label=L"$V^{\rm (d)}_2$")
# # #     # # axs[1,1].plot(r,Objects["potential"][:,2,2],color="red","--")
# # #     # # axs[1,1].plot(r,Objects["potential"][:,3,3],color="blue","--")

# # #     # # inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25])  # [left, bottom, width, height]

# # #     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,2,2], color="green")
# # #     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,3,3], color="orange")
# # #     # # inset_ax.plot(r,Objects["potential"][:,2,2],color="red","--")
# # #     # # inset_ax.plot(r,Objects["potential"][:,3,3],color="blue","--")

# # #     # # inset_ax.set_xlim(1.8, 2.1)  # Zoom into x-range
# # #     # # inset_ax.set_ylim(50000,63000)  # Zoom into y-range

# # #     # # axs[1,1].plot([1.8,1.8],[5e4,6.3e4], color="grey"     ,alpha=0.5)
# # #     # # axs[1,1].plot([2.1,2.1],[5e4,6.3e4], color="grey"     ,alpha=0.5)
# # #     # # axs[1,1].plot([1.8,2.1],[5e4,5e4], color="grey"       ,alpha=0.5)
# # #     # # axs[1,1].plot([1.8,2.1],[6.3e4,6.3e4], color="grey"   ,alpha=0.5)
# # #     # # axs[1,1].plot([2.1,2.668],[5e4,5e4], color="grey"     ,alpha=0.5)
# # #     # # axs[1,1].plot([2.1,2.668],[6.3e4,7.09e4], color="grey",alpha=0.5)

# # #     # # inset_ax.set_xticklabels([])  # Remove x-axis numbers
# # #     # # inset_ax.set_yticklabels([])  # Remove y-axis numbers

# # #     # # # axs[1,1].set_xlabel("Bond Length")
# # #     # # axs[1,1].set_ylabel(L"Potential Energy, $cm^{-1}$")

# # #     # for i=1:dim
# # #     #     for j=i+1:dim
# # #     #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #     #             # axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
# # #     #             axs[2,1].plot(r,Objects["nac"][:,i,j]) #,"--",alpha=0.5)
# # #     #         end
# # #     #     end
# # #     # end

# # #     # # axs[2,1].plot(r,Objects["nac"][:,2,3],"k") #,"--",alpha=0.5)


# # #     # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
# # #     # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$")
# # #     # # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
# # #     # # axs[1,1].legend()
# # #     # # axs[2,1].legend()

# # #     # sf = Calculation["method"].states[end]
# # #     # Emax = Objects["potential"][end,sf,sf]

# # #     # si = Calculation["method"].states[1]
# # #     # Emin = minimum(Objects["potential"][:,si,si])

# # #     # axs[1,1].set_ylim(Emin,1.1*Emax)
# # #     # # axs[1,1].set_xlim(1.35,3.5)

# # #     # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/property_based_pot_guess.png",dpi=300,bbox_inches="tight")
# # # end
# # # plt.figure()

# # # dim=3
# # # for i=1:dim
# # #     for j=i:dim
# # #         plt.plot(r,Diabatic_Objects["potential"][:,i,j])
# # #         plt.plot(r,Objects["potential"][:,i,j],"--")
# # #     end
# # # end

# # # Calculation["method"].states = [2,3]

# # # Potential[2].Lval = r
# # # Potential[2].Rval = Diabatic_Objects["potential"][:,1,1]

# # # include("Build_Hamiltonian_Matrix.jl")


# # # if Calculation["method"].abinitio_fit == true
# # #     fit_abinitio()
# # # elseif Calculation["method"].abinitio_fit == false
# # #     U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
# # #     #
# # #     # fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])

# # #     # plt.subplots_adjust(wspace=0, hspace=0)

# # #     # # axs[1,1].set_title("MRCI aug-cc-pVQZ-X2C | occ = 8330, closed = 5110")
# # #     # for i=1:dim
# # #     #     if lowercase(Calculation["method"].regularisation) == "potential"
# # #     #         if i in Calculation["method"].states
# # #     #             axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
# # #     #             axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
# # #     #         end
# # #     #     else
# # #     #         for j=i:dim
# # #     #             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #     #                 axs[1,1].plot(r,Diabatic_Objects[Calculation["method"].regularisation][:,i,j],label="F"*string(i)*string(j))
# # #     #                 axs[1,1].plot(r,Objects[Calculation["method"].regularisation][:,i,j],"--")
# # #     #             end
# # #     #         end
# # #     #     end
# # #     # end

# # #     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], color="green",label= L"$V^{\rm (d)}_1$")
# # #     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], color="orange",label=L"$V^{\rm (d)}_2$")
# # #     # # axs[1,1].plot(r,Objects["potential"][:,2,2],color="red","--")
# # #     # # axs[1,1].plot(r,Objects["potential"][:,3,3],color="blue","--")

# # #     # # inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25])  # [left, bottom, width, height]

# # #     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,2,2], color="green")
# # #     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,3,3], color="orange")
# # #     # # inset_ax.plot(r,Objects["potential"][:,2,2],color="red","--")
# # #     # # inset_ax.plot(r,Objects["potential"][:,3,3],color="blue","--")

# # #     # # inset_ax.set_xlim(1.8, 2.1)  # Zoom into x-range
# # #     # # inset_ax.set_ylim(50000,63000)  # Zoom into y-range

# # #     # # axs[1,1].plot([1.8,1.8],[5e4,6.3e4], color="grey"     ,alpha=0.5)
# # #     # # axs[1,1].plot([2.1,2.1],[5e4,6.3e4], color="grey"     ,alpha=0.5)
# # #     # # axs[1,1].plot([1.8,2.1],[5e4,5e4], color="grey"       ,alpha=0.5)
# # #     # # axs[1,1].plot([1.8,2.1],[6.3e4,6.3e4], color="grey"   ,alpha=0.5)
# # #     # # axs[1,1].plot([2.1,2.668],[5e4,5e4], color="grey"     ,alpha=0.5)
# # #     # # axs[1,1].plot([2.1,2.668],[6.3e4,7.09e4], color="grey",alpha=0.5)

# # #     # # inset_ax.set_xticklabels([])  # Remove x-axis numbers
# # #     # # inset_ax.set_yticklabels([])  # Remove y-axis numbers

# # #     # # # axs[1,1].set_xlabel("Bond Length")
# # #     # # axs[1,1].set_ylabel(L"Potential Energy, $cm^{-1}$")

# # #     # for i=1:dim
# # #     #     for j=i+1:dim
# # #     #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #     #             # axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
# # #     #             axs[2,1].plot(r,Objects["nac"][:,i,j]) #,"--",alpha=0.5)
# # #     #         end
# # #     #     end
# # #     # end

# # #     # # axs[2,1].plot(r,Objects["nac"][:,2,3],"k") #,"--",alpha=0.5)


# # #     # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
# # #     # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$")
# # #     # # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
# # #     # # axs[1,1].legend()
# # #     # # axs[2,1].legend()

# # #     # sf = Calculation["method"].states[end]
# # #     # Emax = Objects["potential"][end,sf,sf]

# # #     # si = Calculation["method"].states[1]
# # #     # Emin = minimum(Objects["potential"][:,si,si])

# # #     # axs[1,1].set_ylim(Emin,1.1*Emax)
# # #     # # axs[1,1].set_xlim(1.35,3.5)

# # #     # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/property_based_pot_guess.png",dpi=300,bbox_inches="tight")
# # # end

# # # plt.figure()

# # # dim=3
# # # for i=1:dim
# # #     for j=i:dim
# # #         plt.plot(r,Diabatic_Objects["potential"][:,i,j])
# # #         plt.plot(r,Objects["potential"][:,i,j],"--")
# # #     end
# # # end
# # # 

# # # dim = 4
# # # if Calculation["method"].abinitio_fit == true
# # #     fit_abinitio()
# # # elseif Calculation["method"].abinitio_fit == false
# # #     U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
# # #     # #
# # #     # fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])

# # #     # plt.subplots_adjust(wspace=0, hspace=0)

# # #     # # axs[1,1].set_title("MRCI aug-cc-pVQZ-X2C | occ = 8330, closed = 5110")
# # #     # for i=1:dim
# # #     #     if lowercase(Calculation["method"].regularisation) == "potential"
# # #     #         if i in Calculation["method"].states
# # #     #             axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
# # #     #             axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
# # #     #         end
# # #     #     else
# # #     #         for j=i:dim
# # #     #             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #     #                 axs[1,1].plot(r,Diabatic_Objects[Calculation["method"].regularisation][:,i,j],label="F"*string(i)*string(j))
# # #     #                 axs[1,1].plot(r,Objects[Calculation["method"].regularisation][:,i,j],"--")
# # #     #             end
# # #     #         end
# # #     #     end
# # #     # end

# # #     # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], color="green",label= L"$V^{\rm (d)}_1$")
# # #     # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], color="orange",label=L"$V^{\rm (d)}_2$")
# # #     # axs[1,1].plot(r,Objects["potential"][:,2,2],color="red","--")
# # #     # axs[1,1].plot(r,Objects["potential"][:,3,3],color="blue","--")

# # #     # inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25])  # [left, bottom, width, height]

# # #     # inset_ax.plot(r,Diabatic_Objects["potential"][:,2,2], color="green")
# # #     # inset_ax.plot(r,Diabatic_Objects["potential"][:,3,3], color="orange")
# # #     # inset_ax.plot(r,Objects["potential"][:,2,2],color="red","--")
# # #     # inset_ax.plot(r,Objects["potential"][:,3,3],color="blue","--")

# # #     # inset_ax.set_xlim(1.8, 2.1)  # Zoom into x-range
# # #     # inset_ax.set_ylim(50000,63000)  # Zoom into y-range

# # #     # axs[1,1].plot([1.8,1.8],[5e4,6.3e4], color="grey"     ,alpha=0.5)
# # #     # axs[1,1].plot([2.1,2.1],[5e4,6.3e4], color="grey"     ,alpha=0.5)
# # #     # axs[1,1].plot([1.8,2.1],[5e4,5e4], color="grey"       ,alpha=0.5)
# # #     # axs[1,1].plot([1.8,2.1],[6.3e4,6.3e4], color="grey"   ,alpha=0.5)
# # #     # axs[1,1].plot([2.1,2.668],[5e4,5e4], color="grey"     ,alpha=0.5)
# # #     # axs[1,1].plot([2.1,2.668],[6.3e4,7.09e4], color="grey",alpha=0.5)

# # #     # inset_ax.set_xticklabels([])  # Remove x-axis numbers
# # #     # inset_ax.set_yticklabels([])  # Remove y-axis numbers

# # #     # # axs[1,1].set_xlabel("Bond Length")
# # #     # # axs[1,1].set_ylabel(L"Potential Energy, $cm^{-1}$")

# # #     # for i=1:dim
# # #     #     for j=i+1:dim
# # #     #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #     #             # axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
# # #     #             axs[2,1].plot(r,Objects["nac"][:,i,j]) #,"--",alpha=0.5)
# # #     #         end
# # #     #     end
# # #     # end

# # #     # # axs[2,1].plot(r,Objects["nac"][:,2,3],"k") #,"--",alpha=0.5)


# # #     # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
# # #     # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$")
# # #     # # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
# # #     # # axs[1,1].legend()
# # #     # # axs[2,1].legend()

# # #     # sf = Calculation["method"].states[end]
# # #     # Emax = Objects["potential"][end,sf,sf]

# # #     # si = Calculation["method"].states[1]
# # #     # Emin = minimum(Objects["potential"][:,si,si])

# # #     # axs[1,1].set_ylim(Emin,1.1*Emax)
# # #     # axs[1,1].set_xlim(1.35,3.5)

# # #     # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/property_based_pot_guess.png",dpi=300,bbox_inches="tight")
# # # end

# # # Calculation["method"].diabatisation = "forward-evolution"

# # # U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))

# # # fig, axs = plt.subplots(2, 2,sharex=true)
# # # #
# # # axs[1,1].tick_params(axis="both", labelsize=8)
# # # axs[1,2].tick_params(axis="both", labelsize=8)
# # # axs[2,1].tick_params(axis="both", labelsize=8)
# # # axs[2,2].tick_params(axis="both", labelsize=8)
# # # lw = 0.9
# # # #
# # # axs[1,1].set_title("Forward Evolution",fontsize=8)

# # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,1,1], linewidth=lw,color="black",label = L"V$^{\rm(d)}_{1}$")
# # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], linewidth=lw,color="grey",label = L"V$^{\rm(d)}_{2}$")
# # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], linewidth=lw,color="orange",label = L"V$^{\rm(d)}_{3}$")
# # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,4,4], linewidth=lw,color="brown",label = L"V$^{\rm(d)}_{4}$")

# # # axs[1,1].plot(r,Objects["potential"][:,1,1], "--", alpha=0.5,linewidth=lw, color="red",    label = L"$C^{2}\Sigma^{+}+\epsilon K_{11}$")
# # # axs[1,1].plot(r,Objects["potential"][:,2,2], "--", alpha=0.5,linewidth=lw, color="green",   label = L"$(2)^{2}\Sigma^{+}+\epsilon K_{22}$")
# # # axs[1,1].plot(r,Objects["potential"][:,3,3], "--", alpha=0.5,linewidth=lw, color="blue",   label = L"$(3)^{2}\Sigma^{+}+\epsilon K_{33}$")
# # # axs[1,1].plot(r,Objects["potential"][:,4,4], "--", alpha=0.5,linewidth=lw, color="magenta",label = L"$(4)^{2}\Sigma^{+}+\epsilon K_{44}$")

# # # axs[2,1].plot(r,Diabatic_Objects["potential"][:,2,3],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{23}$")
# # # axs[2,1].plot(r,Diabatic_Objects["potential"][:,1,3],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{13}$")
# # # axs[2,1].plot(r,Diabatic_Objects["potential"][:,1,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{14}$")
# # # axs[2,1].plot(r,Diabatic_Objects["potential"][:,2,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{24}$")
# # # axs[2,1].plot(r,Diabatic_Objects["potential"][:,3,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{34}$")
# # # axs[2,1].plot(r,Diabatic_Objects["potential"][:,1,2],linewidth=lw,color="red",label=L"W$^{\rm(1)}_{12}$")



# # # Calculation["method"].diabatisation = "backward-evolution"

# # # U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))


# # # axs[1,2].set_title("Backward Evolution",fontsize=8)

# # # axs[1,2].plot(r,Diabatic_Objects["potential"][:,1,1], linewidth=lw,color="black",label = L"V$^{\rm(d)}_{1}$")
# # # axs[1,2].plot(r,Diabatic_Objects["potential"][:,2,2], linewidth=lw,color="grey",label = L"V$^{\rm(d)}_{2}$")
# # # axs[1,2].plot(r,Diabatic_Objects["potential"][:,3,3], linewidth=lw,color="orange",label = L"V$^{\rm(d)}_{3}$")
# # # axs[1,2].plot(r,Diabatic_Objects["potential"][:,4,4], linewidth=lw,color="brown",label = L"V$^{\rm(d)}_{4}$")

# # # axs[1,2].plot(r,Objects["potential"][:,1,1], "--", linewidth=lw, alpha=0.5,color="red",    label = L"$C^{2}\Sigma^{+}+\epsilon K_{11}$")
# # # axs[1,2].plot(r,Objects["potential"][:,2,2], "--", linewidth=lw, alpha=0.5,color="green",   label = L"$(2)^{2}\Sigma^{+}+\epsilon K_{22}$")
# # # axs[1,2].plot(r,Objects["potential"][:,3,3], "--", linewidth=lw, alpha=0.5,color="blue",   label = L"$(3)^{2}\Sigma^{+}+\epsilon K_{33}$")
# # # axs[1,2].plot(r,Objects["potential"][:,4,4], "--", linewidth=lw, alpha=0.5,color="magenta",label = L"$(4)^{2}\Sigma^{+}+\epsilon K_{44}$")

# # # axs[2,2].plot(r,Diabatic_Objects["potential"][:,2,3],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{23}$")
# # # axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,3],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{13}$")
# # # axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{14}$")
# # # axs[2,2].plot(r,Diabatic_Objects["potential"][:,2,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{24}$")
# # # axs[2,2].plot(r,Diabatic_Objects["potential"][:,3,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{34}$")
# # # axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,2],linewidth=lw,color="red",label=L"W$^{\rm(1)}_{12}$")


# # # #
# # # axs[2,1].set_xlabel(L"Bond Length R, $\mathrm{\AA}$",fontsize=8)
# # # axs[2,2].set_xlabel(L"Bond Length R, $\mathrm{\AA}$",fontsize=8)
# # # #
# # # axs[2,1].set_ylabel(L"DC, $\rm cm^{-1}$",fontsize=8)
# # # #
# # # axs[1,1].set_ylabel(L"Potential, cm$^{-1}$",fontsize=8)
# # # #
# # # axs[1,1].set_ylim(25000,120000)
# # # axs[1,2].set_ylim(25000,120000)
# # # axs[2,1].set_ylim(-10100,21000)
# # # axs[2,2].set_ylim(-10100,21000)
# # # axs[1,1].set_xlim(r[1],6)
# # # axs[2,1].set_xlim(r[1],6)
# # # axs[1,2].set_xlim(r[1],6)
# # # axs[2,2].set_xlim(r[1],6)

# # # axs[1,2].set_yticklabels([])
# # # axs[2,2].set_yticklabels([])

# # # plt.subplots_adjust(hspace=0.0,wspace=0.1)
# # # fig.set_size_inches(6.69, 3.1)

# # # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/N-level/Evolution/CH/plots/paper/CH_forward_backward_evo_diff_PECs_DCs_simplified.png",bbox_inches="tight",dpi=600)









# # # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/CH_scale_NAC.png",bbox_inches="tight",dpi=600)
# # # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/YO_shift_NAC_2.png",bbox_inches="tight",dpi=600)

# # # plt.figure()

# # # for i=1:dim
# # #     for j=i:dim
# # #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #             plt.plot(r,Objects["nac"][:,i,j])
# # #             # plt.plot(r,[Objects["regularised_nac"][idx][i,j] for idx=1:lastindex(r)],"--")
# # #         end
# # #     end
# # # end

# # # save_diabatisation(Objects, Diabatic_Objects, lowercase(Calculation["method"].diabatisation), input_properties, fname, special_name = "4paper")



# # # fig, axs = plt.subplots(2, 2,sharex=true)
# # # axs[1,1].tick_params(axis="both", labelsize=8)
# # # axs[1,2].tick_params(axis="both", labelsize=8)
# # # axs[2,1].tick_params(axis="both", labelsize=8)
# # # axs[2,2].tick_params(axis="both", labelsize=8)
# # # lw = 0.9
# # # #
# # # ## diabats
# # # axs[1,2].text(5.2,110000, L"\rm CH",fontsize=10,color="black")

# # # axs[1,2].plot(r,Diabatic_Objects["potential"][:,1,1], linewidth=lw, color="black",label = L"V$^{\rm(d)}_{1}$")
# # # axs[1,2].plot(r,Diabatic_Objects["potential"][:,2,2], linewidth=lw, color="grey",label = L"V$^{\rm(d)}_{2}$")
# # # axs[1,2].plot(r,Diabatic_Objects["potential"][:,3,3], linewidth=lw, color="orange",label = L"V$^{\rm(d)}_{3}$")
# # # axs[1,2].plot(r,Diabatic_Objects["potential"][:,4,4], linewidth=lw, color="blue",label = L"V$^{\rm(d)}_{4}$")


# # # axs[1,2].text(1.83,59100,  L"V$^{\rm(d)}_{3}$",fontsize=8,color="orange")
# # # axs[1,2].text(1.58 ,31000,  L"V$^{\rm(d)}_{1}$",fontsize=8,color="black")
# # # axs[1,2].text(3.15,108000,  L"V$^{\rm(d)}_{2}$",fontsize=8,color="grey")
# # # axs[1,2].text(2.87,86300,  L"V$^{\rm(d)}_{4}$",fontsize=8,color="blue")

# # # axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,2], linewidth=lw, color="blue", label=L"$\mathcal{D}_{12}$")
# # # axs[2,2].plot(r,Diabatic_Objects["potential"][:,2,3], linewidth=lw, color="red",  label=L"$\mathcal{D}_{23}$")
# # # axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,3], linewidth=lw, color="green",label=L"$\mathcal{D}_{13}$")
# # # axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,4], linewidth=lw, color="black",label=L"$\mathcal{D}_{14}$")
# # # axs[2,2].plot(r,Diabatic_Objects["potential"][:,2,4], linewidth=lw, color="orange",label=L"$\mathcal{D}_{24}$")
# # # axs[2,2].plot(r,Diabatic_Objects["potential"][:,3,4], linewidth=lw, color="purple",label=L"$\mathcal{D}_{34}$")


# # # axs[2,2].text(1.25,3000,  L"$V^{\rm(d)}_{12}$",fontsize=8,color="blue")
# # # axs[2,2].text(2.33, 10000, L"$V^{\rm(d)}_{23}$",fontsize=8,color="red")
# # # axs[2,2].text(3,-4500 , L"$V^{\rm(d)}_{13}$",fontsize=8,color="green")
# # # axs[2,2].text(4.5,3000, L"$V^{\rm(d)}_{14}$",fontsize=8,color="black")
# # # axs[2,2].text(3.2 ,4300,  L"$V^{\rm(d)}_{24}$",fontsize=8,color="orange")
# # # axs[2,2].text(2.25,-2500,  L"$V^{\rm(d)}_{34}$",fontsize=8,color="purple")

# # # #
# # # ## adiabats

# # # atom1, atom2 = Calculation["method"].atoms
# # # kinetic_factor =  KE_factor(atom1, atom2)


# # # axs[1,1].plot(r,Objects["potential"][:,1,1] .+ (kinetic_factor .* [k[1,1] for k in Objects["K_matrix"]]), linewidth=lw, color="red",    label = L"$C^{2}\Sigma^{+}+\epsilon K_{11}$")
# # # axs[1,1].plot(r,Objects["potential"][:,2,2] .+ (kinetic_factor .* [k[2,2] for k in Objects["K_matrix"]]), linewidth=lw, color="green",   label = L"$(2)^{2}\Sigma^{+}+\epsilon K_{22}$")
# # # axs[1,1].plot(r,Objects["potential"][:,3,3] .+ (kinetic_factor .* [k[3,3] for k in Objects["K_matrix"]]), linewidth=lw, color="blue",   label = L"$(3)^{2}\Sigma^{+}+\epsilon K_{33}$")
# # # axs[1,1].plot(r,Objects["potential"][:,4,4] .+ (kinetic_factor .* [k[4,4] for k in Objects["K_matrix"]]), linewidth=lw, color="magenta",label = L"$(4)^{2}\Sigma^{+}+\epsilon K_{44}$")


# # # axs[1,1].text(3.5,33000,  L"$V_{C^{2}\Sigma^{+}}+\epsilon K_{11}$",fontsize=8,color="red")
# # # axs[1,1].text(3.5,44000,  L"$V_{2^{2}\Sigma^{+}}+\epsilon K_{22}$",fontsize=8,color="green")
# # # axs[1,1].text(3.5,72000,  L"$V_{3^{2}\Sigma^{+}}+\epsilon K_{33}$",fontsize=8,color="blue")
# # # axs[1,1].text(3.5,110000, L"$V_{4^{2}\Sigma^{+}}+\epsilon K_{44}$",fontsize=8,color="magenta")


# # # axs[2,1].plot(r,[w[1,2] for w in Objects["regularised_nac"]], linewidth=lw,color="blue",label=L"W$^{\rm(1)}_{12}$")
# # # axs[2,1].plot(r,[w[2,3] for w in Objects["regularised_nac"]], linewidth=lw,color="red",label=L"W$^{\rm(1)}_{23}$")
# # # axs[2,1].plot(r,[w[1,3] for w in Objects["regularised_nac"]], linewidth=lw,color="green",label=L"W$^{\rm(1)}_{13}$")
# # # axs[2,1].plot(r,[w[1,4] for w in Objects["regularised_nac"]], linewidth=lw,color="black",label=L"W$^{\rm(1)}_{14}$")
# # # axs[2,1].plot(r,[w[2,4] for w in Objects["regularised_nac"]], linewidth=lw,color="orange",label=L"W$^{\rm(1)}_{24}$")
# # # axs[2,1].plot(r,[w[3,4] for w in Objects["regularised_nac"]], linewidth=lw,color="purple",label=L"W$^{\rm(1)}_{34}$")

# # # axs[2,1].text(2.95, 1.4  , L"$\langle C^{2}\Sigma^{+}|\frac{d}{dr}|2^{2}\Sigma^{+}\rangle$",fontsize=8,color="blue")
# # # axs[2,1].text(2.95, 1.05  ,L"$\langle 2^{2}\Sigma^{+}|\frac{d}{dr}|3^{2}\Sigma^{+}\rangle$",fontsize=8,color="red")
# # # axs[2,1].text(2.95, 0.35 , L"$\langle C^{2}\Sigma^{+}|\frac{d}{dr}|3^{2}\Sigma^{+}\rangle$",fontsize=8,color="green")
# # # axs[2,1].text(2.95, -0.25 , L"$\langle C^{2}\Sigma^{+}|\frac{d}{dr}|4^{2}\Sigma^{+}\rangle$",fontsize=8,color="black")
# # # axs[2,1].text(2.95, 0.7  , L"$\langle 2^{2}\Sigma^{+}|\frac{d}{dr}|4^{2}\Sigma^{+}\rangle$",fontsize=8,color="orange")
# # # axs[2,1].text(2.95, 1.75 , L"$\langle 3^{2}\Sigma^{+}|\frac{d}{dr}|4^{2}\Sigma^{+}\rangle$",fontsize=8,color="purple")
# # # #
# # # axs[2,1].set_xlabel(L"Bond Length R, $\mathrm{\AA}$",fontsize=8)
# # # axs[2,2].set_xlabel(L"Bond Length R, $\mathrm{\AA}$",fontsize=8)
# # # #
# # # axs[2,1].set_ylabel("Coupling",fontsize=8)
# # # #
# # # axs[1,1].set_ylabel(L"Potential Energy, cm$^{-1}$",fontsize=8)
# # # #
# # # axs[1,1].set_ylim(25000,120000)
# # # axs[1,2].set_ylim(25000,120000)
# # # axs[1,1].set_xlim(r[1],6)
# # # axs[2,1].set_xlim(r[1],6)

# # # plt.subplots_adjust(hspace=0.0,wspace=0.4)
# # # fig.set_size_inches(4.5,5)

# # # axs[1,1].tick_params(axis="both", labelsize=8)
# # # axs[2,1].tick_params(axis="both", labelsize=8)
# # # axs[1,2].tick_params(axis="both", labelsize=8)
# # # axs[2,2].tick_params(axis="both", labelsize=8)

# # # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/wavefunction_comparer/CH/plots/CH_diabatisation_adi-dia_4x5.png",bbox_inches="tight",dpi=600)





# # #
# # ######## CO

# # # Va2 = deepcopy(Potential[2].Rval)
# # # rVa2 = deepcopy(Potential[2].Lval)


# # # QD = deepcopy(Diabatic_Objects["potential"][:,3,3])
# # # R = deepcopy(r)

# # # Calculation["method"].states = [1,2]

# # # Potential[2].Lval = r
# # # Potential[2].Rval = QD

# # # if Calculation["method"].abinitio_fit == true
# # #     fit_abinitio()
# # # elseif Calculation["method"].abinitio_fit == false
# # #     U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
# # #     #
# # #     # fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])

# # #     # plt.subplots_adjust(wspace=0, hspace=0)

# # #     # # axs[1,1].set_title("MRCI aug-cc-pVQZ-X2C | occ = 8330, closed = 5110")
# # #     # for i=1:dim
# # #     #     if lowercase(Calculation["method"].regularisation) == "potential"
# # #     #         if i in Calculation["method"].states
# # #     #             axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
# # #     #             axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
# # #     #         end
# # #     #     else
# # #     #         for j=i:dim
# # #     #             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #     #                 axs[1,1].plot(r,Diabatic_Objects[Calculation["method"].regularisation][:,i,j],label="F"*string(i)*string(j))
# # #     #                 axs[1,1].plot(r,Objects[Calculation["method"].regularisation][:,i,j],"--")
# # #     #             end
# # #     #         end
# # #     #     end
# # #     # end

# # #     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], color="green",label= L"$V^{\rm (d)}_1$")
# # #     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], color="orange",label=L"$V^{\rm (d)}_2$")
# # #     # # axs[1,1].plot(r,Objects["potential"][:,2,2],color="red","--")
# # #     # # axs[1,1].plot(r,Objects["potential"][:,3,3],color="blue","--")

# # #     # # inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25])  # [left, bottom, width, height]

# # #     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,2,2], color="green")
# # #     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,3,3], color="orange")
# # #     # # inset_ax.plot(r,Objects["potential"][:,2,2],color="red","--")
# # #     # # inset_ax.plot(r,Objects["potential"][:,3,3],color="blue","--")

# # #     # # inset_ax.set_xlim(1.8, 2.1)  # Zoom into x-range
# # #     # # inset_ax.set_ylim(50000,63000)  # Zoom into y-range

# # #     # # axs[1,1].plot([1.8,1.8],[5e4,6.3e4], color="grey"     ,alpha=0.5)
# # #     # # axs[1,1].plot([2.1,2.1],[5e4,6.3e4], color="grey"     ,alpha=0.5)
# # #     # # axs[1,1].plot([1.8,2.1],[5e4,5e4], color="grey"       ,alpha=0.5)
# # #     # # axs[1,1].plot([1.8,2.1],[6.3e4,6.3e4], color="grey"   ,alpha=0.5)
# # #     # # axs[1,1].plot([2.1,2.668],[5e4,5e4], color="grey"     ,alpha=0.5)
# # #     # # axs[1,1].plot([2.1,2.668],[6.3e4,7.09e4], color="grey",alpha=0.5)

# # #     # # inset_ax.set_xticklabels([])  # Remove x-axis numbers
# # #     # # inset_ax.set_yticklabels([])  # Remove y-axis numbers

# # #     # # # axs[1,1].set_xlabel("Bond Length")
# # #     # # axs[1,1].set_ylabel(L"Potential Energy, $cm^{-1}$")

# # #     # for i=1:dim
# # #     #     for j=i+1:dim
# # #     #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
# # #     #             # axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
# # #     #             axs[2,1].plot(r,Objects["nac"][:,i,j]) #,"--",alpha=0.5)
# # #     #         end
# # #     #     end
# # #     # end

# # #     # # axs[2,1].plot(r,Objects["nac"][:,2,3],"k") #,"--",alpha=0.5)


# # #     # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
# # #     # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$")
# # #     # # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
# # #     # # axs[1,1].legend()
# # #     # # axs[2,1].legend()

# # #     # sf = Calculation["method"].states[end]
# # #     # Emax = Objects["potential"][end,sf,sf]

# # #     # si = Calculation["method"].states[1]
# # #     # Emin = minimum(Objects["potential"][:,si,si])

# # #     # axs[1,1].set_ylim(Emin,1.1*Emax)
# # #     # # axs[1,1].set_xlim(1.35,3.5)

# # #     # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/property_based_pot_guess.png",dpi=300,bbox_inches="tight")
# # # end

# # # rep = deepcopy(Diabatic_Objects["potential"][:,2,2])



# # # f = NonAdiabaticCoupling[[1,2]].type
# # # p =  NonAdiabaticCoupling[[1,2]].Rval
# # # #
# # # ##
# # # println("CATS")
# # # mixing_angle_ij = compute_mixing_angle(r, f, [1,2], p)
# # # #
# # # c = cos.(mixing_angle_ij)
# # # s = sin.(mixing_angle_ij)
# # # #
# # # U = zeros(Float64,length(r),2,2)
# # # U[:,1,1] .= c
# # # U[:,2,2] .= c
# # # U[:,1,2] .= -s
# # # U[:,2,1] .= s

# # # Va = zeros(Float64,length(r),2,2)
# # # Va[:,1,1] .= Objects["potential"][:,1,1]
# # # Va[:,2,2] .= QD

# # # Vd = map(idx -> U[idx,:,:]' * Va[idx,:,:] * U[idx,:,:], collect(1:lastindex(r)))

# # # rep = [v[2,2] for v in Vd]


# # # Calculation["grid"].npoints = 3001
# # # include("Build_Hamiltonian_Matrix.jl")

# # # fig, axs = plt.subplots(2,2,figsize=[6,4],sharex=true)
# # # plt.subplots_adjust(hspace=0, wspace=0.65)  # Remove vertical space


# # # # plt.plot(Potential[5].Lval,Potential[5].Rval)
# # # # plt.plot(Potential[6].Lval,Potential[6].Rval)
# # # # plt.plot(Potential[7].Lval,Potential[7].Rval)

# # # axs[1,1].plot(Potential[1].Lval,Potential[1].Rval,"blue")
# # # axs[1,1].plot(Potential[3].Lval,Potential[3].Rval,"orange")
# # # axs[1,1].plot(rVa2,Va2,"limegreen")
# # # # axs[1,1].plot(Potential[4].Lval,Potential[4].Rval)
# # # axs[1,1].plot(R,QD,color="red","--")
# # # axs[1,1].set_ylabel(L"Potential Energy, cm$^{-1}$")
# # # axs[1,1].text(1.4,92000,L"$B^1\Sigma^+$",fontsize=9,color="blue")
# # # axs[1,1].text(1.06,97000,L"$C^1\Sigma^+$",fontsize=9,color="limegreen")
# # # axs[1,1].text(1.1,113000,L"$(IV)^1\Sigma^+$",fontsize=9,color="orange")
# # # axs[1,1].text(1.4,100000,L"$QD^1\Sigma^+$",fontsize=9,color="red")


# # # axs[1,1].set_ylim(85000,120000)


# # # global r
# # # r = LinRange(Calculation["grid"].range[1],
# # #              Calculation["grid"].range[2],
# # #              Calculation["grid"].npoints)  

# # # axs[2,1].plot(r,Objects["nac"][:,1,2],color="k")
# # # axs[2,1].plot(r,Objects["nac"][:,2,3],color="red")
# # # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$" )
# # # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
# # # axs[2,1].text(0.95,45,L"$\langle C^1\Sigma^+| d/dr | (IV)^1\Sigma^+ \rangle$",color="red",fontsize=9)
# # # axs[2,1].text(1.1,40,L"$\langle B^1\Sigma^+| d/dr | QD^1\Sigma^+ \rangle$",color="k",fontsize=9)
# # # axs[2,1].set_ylim(0,50)

# # # axs[1,2].plot(R,rep,"orange",alpha=0.25)
# # # axs[1,2].plot(Potential[5].Lval,Potential[5].Rval,"blue")
# # # axs[1,2].plot(Potential[6].Lval,Potential[6].Rval,"green")
# # # axs[1,2].set_ylabel(L"Potential Energy, cm$^{-1}$")
# # # axs[1,2].text(1.175,105000,L"$C^1\Sigma^{+\rm(d)}$",fontsize=9,color="green")
# # # axs[1,2].text(1.4,95000,L"$B^1\Sigma^{+\rm(d)}$",fontsize=9,color="blue")
# # # axs[1,2].text(1.1,115000,L"$(3)^1\Sigma^{+\rm(d)}$",fontsize=9,color="orange")



# # # axs[2,2].plot(Potential[8].Lval,Potential[8].Rval.*219474.6313708000,color="black")
# # # axs[2,2].plot(Potential[9].Lval,Potential[9].Rval,color="red")
# # # axs[2,2].set_ylabel(L"DC, cm$^{-1}$")
# # # axs[2,2].set_xlabel(L"Bond Length, $\rm \AA$")
# # # axs[2,2].set_ylim(0,4000)

# # # axs[2,2].text(0.95,3550,L"$\langle C^1\Sigma^{+\rm(d)}| V^{\rm(d)} | QD^1\Sigma^+ \rangle$",color="red",fontsize=9)
# # # axs[2,2].text(0.95,3000,L"$\langle B^1\Sigma^{+\rm(d)}| V^{\rm(d)} | (3)^1\Sigma^{+ \rm (d)} \rangle$",color="k",fontsize=9)
# # # # axs[2,2].text(1.4,500,L"$\times 50000$",color="k",fontsize=8)


# # # axs[1,2].set_xlim(0.9,1.6)
# # # axs[1,2].set_ylim(85000,120000)

# # # plt.savefig("./CO_2state_diabatisation_2.png", bbox_inches="tight",dpi=300)







# # # function mixang_accumulation_metric(gamma,L)
# # #     x = L/gamma
# # #     return log(x+sqrt(1+x^2))/x
# # # end

# # # g = LinRange(0.001,6,2000)

# # # Ls = [1,2,3,4,5]

# # # plt.figure()

# # # for L in Ls
# # #     plt.plot(g,mixang_accumulation_metric.(g,L),label="FC width = $L")
# # # end
# # # plt.legend()
# # # plt.xlabel(L"NAC width, $\gamma$ ($\AA$)")
# # # plt.ylabel("Accumulated Mixing Factor")
# # # plt.xscale("log")

# # # plt.savefig("AMF.png",dpi=300,bbox_inches="tight")