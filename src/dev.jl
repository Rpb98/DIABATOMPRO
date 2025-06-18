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
read_file("../Supplementary/N2/DIABATOM-PRO_N2.inp")
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
# if Calculation["method"].abinitio_fit == true
#     fit_abinitio()
# end

#
## compute vibronic energies and wavefunctions for a non-rotating molecule if
## the 'vibronic_solver' key is in Calculation
# if haskey(Calculation, "vibronic_solver")
#     #
#     @time contr_vib_wfn, E_vib_contr = vibronic_eigensolver(collect(r), PotMat, Calculation["grid"].npoints, collect(keys(Potential)), Calculation["method"].atoms...)
# end



function serialise_electronic_vib_indices(state, v_idx, contraction_array)
    return ( state - 1 ) * contraction_array[ state ] + v_idx
end
#
vmax = Calculation["vibronic_solver"].contraction
contracted_vibronic_dim = sum(vmax)
#
T = zeros(Float64, contracted_vibronic_dim, contracted_vibronic_dim)
V = zeros(Float64, contracted_vibronic_dim, contracted_vibronic_dim)

#
states = collect(keys(Potential))
#
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
H = T + V

eig = eigen(H)


# U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))

# plt.figure()
# for key in keys(Dipole)
#     i, j = key
#     #
#     plt.plot(r,Diabatic_Objects["dipole"][:,i,j])
#     #
#     plt.plot(r,Objects["dipole"][:,i,j],"--")
# end
    #


# plt.ylim(40000,60000)

# plt.figure()
# plt.plot(r, PotMat[:,3,3])
# plt.plot(abinitio[("poten",3)].Lval, abinitio[("poten",3)].Rval,"--")
# # plt.plot(r, PotMat[:,4,4])

# # plt.plot(r, PotMat[:,1,1],"--")
# # plt.plot(r, PotMat[:,2,2],"--")

# # plt.plot(r, PotMat[:,5,5])
# # plt.plot(r, PotMat[:,6,6])
# # plt.plot(r, PotMat[:,7,7])

# plt.ylim(35000,75000)






# #
# ## initialise diabatom object
# diabatom = Dict()
# diabatom["r"] = r
# diabatom["Hamiltonian"] = Hamiltonian
# diabatom["Objects"] = Objects
# diabatom["Calculation"] = Calculation
# diabatom["Potential"] = Potential
# diabatom["SpinOrbit"] = SpinOrbit
# diabatom["EAMC"] = EAMC
# diabatom["Dipole"] = Dipole
# diabatom["NonAdiabaticCoupling"] = NonAdiabaticCoupling

# # run the diabatiser
# if Calculation["method"].abinitio_fit == true
#     fit_abinitio()
# elseif Calculation["method"].abinitio_fit == false
#     U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
    
#     # fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])

#     # plt.subplots_adjust(wspace=0, hspace=0)

#     # # axs[1,1].set_title("MRCI aug-cc-pVQZ-X2C | occ = 8330, closed = 5110")
#     # for i=1:dim
#     #     if lowercase(Calculation["method"].regularisation) == "potential"
#     #         if i in Calculation["method"].states
#     #             axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
#     #             axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
#     #         end
#     #     else
#     #         for j=i:dim
#     #             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#     #                 axs[1,1].plot(r,Diabatic_Objects[Calculation["method"].regularisation][:,i,j],label="F"*string(i)*string(j))
#     #                 axs[1,1].plot(r,Objects[Calculation["method"].regularisation][:,i,j],"--")
#     #             end
#     #         end
#     #     end
#     # end

#     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], color="green",label= L"$V^{\rm (d)}_1$")
#     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], color="orange",label=L"$V^{\rm (d)}_2$")
#     # # axs[1,1].plot(r,Objects["potential"][:,2,2],color="red","--")
#     # # axs[1,1].plot(r,Objects["potential"][:,3,3],color="blue","--")

#     # # inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25])  # [left, bottom, width, height]

#     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,2,2], color="green")
#     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,3,3], color="orange")
#     # # inset_ax.plot(r,Objects["potential"][:,2,2],color="red","--")
#     # # inset_ax.plot(r,Objects["potential"][:,3,3],color="blue","--")

#     # # inset_ax.set_xlim(1.8, 2.1)  # Zoom into x-range
#     # # inset_ax.set_ylim(50000,63000)  # Zoom into y-range

#     # # axs[1,1].plot([1.8,1.8],[5e4,6.3e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.1],[5e4,6.3e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([1.8,2.1],[5e4,5e4], color="grey"       ,alpha=0.5)
#     # # axs[1,1].plot([1.8,2.1],[6.3e4,6.3e4], color="grey"   ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.668],[5e4,5e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.668],[6.3e4,7.09e4], color="grey",alpha=0.5)

#     # # inset_ax.set_xticklabels([])  # Remove x-axis numbers
#     # # inset_ax.set_yticklabels([])  # Remove y-axis numbers

#     # # # axs[1,1].set_xlabel("Bond Length")
#     # # axs[1,1].set_ylabel(L"Potential Energy, $cm^{-1}$")

#     # for i=1:dim
#     #     for j=i+1:dim
#     #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#     #             # axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
#     #             axs[2,1].plot(r,Objects["nac"][:,i,j]) #,"--",alpha=0.5)
#     #         end
#     #     end
#     # end

#     # # axs[2,1].plot(r,Objects["nac"][:,2,3],"k") #,"--",alpha=0.5)


#     # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
#     # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$")
#     # # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
#     # # axs[1,1].legend()
#     # # axs[2,1].legend()

#     # sf = Calculation["method"].states[end]
#     # Emax = Objects["potential"][end,sf,sf]

#     # si = Calculation["method"].states[1]
#     # Emin = minimum(Objects["potential"][:,si,si])

#     # axs[1,1].set_ylim(Emin,1.1*Emax)
#     # axs[1,1].set_xlim(1.35,3.5)

#     # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/property_based_pot_guess.png",dpi=300,bbox_inches="tight")
# end

# fig, axs = plt.subplots(2,2,sharex=true,figsize=[5,3])
# plt.subplots_adjust(wspace=0.4, hspace=0)
# for i=1:dim
#     axs[1,1].plot(r,Objects["potential"][:,i,i])
#     #
#     axs[1,2].plot(r,Diabatic_Objects["potential"][:,i,i])
#     #
#     for j=i+1:dim
#         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#             # axs[2,1].plot(r,Objects["nac"][:,i,j], alpha=0.5)
#             axs[2,1].plot(r,[Objects["regularised_nac"][idx][i,j] for idx=1:lastindex(r)])
#             axs[2,2].plot(r,Diabatic_Objects["potential"][:,i,j])
#         end
#     end
# end

# col = ["k","red","orange","salmon","brown","blue","green"]

# col2 = ["blue","green","orange"]
# fig, axs = plt.subplots(2,1,sharex=true,figsize=[5,4])
# plt.subplots_adjust(wspace=0.4, hspace=0)

# global count,count2
# count = 0
# count2 = 0
# for i=1:dim
#     global count, count2
#     count+=1
#     axs[1,1].plot(r,Objects["potential"][:,i,i],linewidth=1,color=col[count],label=Potential[i].name)
#     #
#     for j=i+1:dim
#         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#             global count2
#             count2+=1
#             axs[2,1].plot(r,Objects["nac"][:,i,j],color=col2[count2],label=string(i)*string(j))

#             Wij_ai = abinitio[("NAC",[i,j])]

#             idx = closest_value_index(Wij_ai.Lval, 2.8)

#             axs[2,1].plot(Wij_ai.Lval[idx:end],Wij_ai.Rval[idx:end],"--",color="k",alpha=1)
#         end
#     end
# end
# # axs[1,1].legend()
# axs[1,1].set_ylim(0,45000)
# axs[1,1].set_xlim(r[1],16)
# axs[1,1].set_ylabel(L"Potential Energy, cm$^{-1}$")
# axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")


# axs[1,1].text(3.85,6500,L"$X^1\Sigma^+$",fontsize=8)
# axs[1,1].text(2.5,11000,L"$(1)^3\Sigma^+$",color="salmon",fontsize=8)
# axs[1,1].text(5.2,18000,L"$(2)^1\Sigma^+$",color="red",fontsize=8)
# axs[1,1].text(4.2,22500,L"$(1)^3\Pi$",color="green",fontsize=8)
# axs[1,1].text(2.6,22000,L"$(1)^1\Pi$",color="blue",fontsize=8)
# axs[1,1].text(3.5,33500,L"$(2)^3\Sigma^+$",color="brown",fontsize=8)
# axs[1,1].text(6.5,33500,L"$(3)^1\Sigma^+$",color="orange",fontsize=8)

# axs[1,1].text(13,15500,"K(4s)+H")
# axs[1,1].text(13,27500,"K(4p)+H")
# axs[1,1].text(13,36500,L"K$^+$+H$^-$",rotation=3)

# # axs[2,1].legend()

# axs[2,1].text(6,-0.2,L"$\langle X^1\Sigma^+|\frac{d}{dr}|(2)^1\Sigma^+\rangle$",color="blue",fontsize=9)
# axs[2,1].text(9.6,0.2,L"$\langle (2)^1\Sigma^+|\frac{d}{dr}|(3)^1\Sigma^+\rangle$",color="green",fontsize=9)
# axs[2,1].text(1.4,0.1,L"$\langle X^1\Sigma^+|\frac{d}{dr}|(3)^1\Sigma^+\rangle$",color="orange",fontsize=9)

# # plt.savefig("./KH_abinitio_PECs_NACs.png",bbox_inches="tight",dpi=300)




# lab_dict = Dict{Int, String}(
#       1 => L"$X^1\Sigma^+$",
#       2 => L"$(2)^1\Sigma^+$",
#       3 => L"$(3)^1\Sigma^+$",
#       4 => L"$(1)^3\Sigma^+$",
#       5 => L"$(2)^3\Sigma^+$",
#       6 => L"$(1)^1\Pi$",
#       7 => L"$(1)^3\Pi$"
#   )
# plt.figure(figsize=[5,3])
# for key in keys(SpinOrbit)
#     i, j = key
#     #
#     if i==j
#         sym = L"$_{\rm z}$"
#     else
#         sym = L"$_{\rm x}$"
#     end
#     #
#     plt.plot(r,Objects["spin-orbit"][:,i,j],label=L"\langle"*lab_dict[i]*L"|SO"*sym*"|"*lab_dict[j]*L"$\rangle$")
# end
# plt.legend(loc="center right",fontsize=8)

# plt.xlabel(L"Bond Length, $\rm \AA$")
# plt.ylabel(L"Spin-Orbit Coupling, cm$^{-1}$")

# plt.savefig("../Supplementary/KH/KH_abinitio_SOCs.png",bbox_inches="tight",dpi=300)


# plt.figure()
# for key in keys(EAMC)
#     i, j = key
#     #
#     plt.plot(r,Objects["lx"][:,i,j])
# end


# plt.figure(figsize=[5,3])
# for key in keys(EAMC)
#     i, j = key
#     #
#     sym = L"$_{\rm x}$"
#     #
#     plt.plot(r,Objects["lx"][:,i,j],label=L"\langle"*lab_dict[i]*L"|L"*sym*"|"*lab_dict[j]*L"$\rangle$")
# end
# plt.legend(loc="center right",fontsize=8)

# plt.xlabel(L"Bond Length, $\rm \AA$")
# plt.ylabel(L"Electronic Angular Momentum, $\hbar$")

# plt.savefig("../Supplementary/KH/KH_abinitio_LX.png",bbox_inches="tight",dpi=300)






# fig, axs = plt.subplots(2,1,figsize=[5,6],gridspec_kw=Dict("height_ratios" => [1.5, 1]))
# plt.subplots_adjust(hspace=0)
# for key in keys(Dipole)
#     i, j = key
#     #
#     name = Dipole[key].name
#     #
#     if occursin("dmx",lowercase(name))
#         C = L"$\mu_{\rm x}$"
#         axs[2,1].plot(r,Objects["dipole"][:,i,j],label=L"\langle"*lab_dict[i]*L"|"*C*"|"*lab_dict[j]*L"$\rangle$")

#     elseif occursin("dmz",lowercase(name))
#         C = L"$\mu_{\rm z}$"
#         axs[1,1].plot(r,Objects["dipole"][:,i,j],label=L"\langle"*lab_dict[i]*L"|"*C*"|"*lab_dict[j]*L"$\rangle$")

#     end
#     #
# end
# axs[1,1].legend(loc="center right",fontsize=8)
# axs[2,1].legend(loc="center right",fontsize=8)

# axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
# axs[2,1].set_ylabel("Dipole Moment, Debye")
# axs[1,1].set_ylabel("Dipole Moment, Debye")


# plt.savefig("../Supplementary/KH/KH_abinitio_DMCs.png",bbox_inches="tight",dpi=300)





# col2 = ["blue","green","orange"]
# plt.figure()
# global count,count2
# count2 = 0
# for i=1:dim
#     global count2
#     #
#     for j=i+1:dim
#         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#             global count2
#             count2+=1
#             #
#             Wij_ai = abinitio[("NAC",[i,j])]
#             idx = closest_value_index(Wij_ai.Lval, 3.5)
#             idxf = closest_value_index(Wij_ai.Lval, 9)

#             #
#             spl = Spline1D(r,Objects["nac"][:,i,j])
#             W_func = spl(Wij_ai.Lval)
#             #
#             # plt.plot(r,Objects["nac"][:,i,j],color=col2[count2],label=string(i)*string(j))

#             diff = (Wij_ai.Rval[idx:idxf] .- W_func[idx:idxf])
#             #
#             RMS = sqrt(sum(diff.^2)/length(diff))
#             println("RMS for ",i," ",j," = ",RMS)
#             plt.plot(Wij_ai.Lval[idx:idxf],(Wij_ai.Rval[idx:idxf] .- W_func[idx:idxf]),color=col2[count2],alpha=1)
#         end
#     end
# end






# fig,axs = plt.subplots(2,2,sharex=true)
# axs[1,1].plot(r,Objects["potential"][:,2,2])
# axs[1,1].plot(r,Objects["potential"][:,3,3])
# axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2],color="k",alpha=0.25)
# axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3],color="k",alpha=0.25)
# # plt.ylim(0.9*minimum(Objects["potential"][:,2,2]),1.1*Objects["potential"][end,3,3])
# axs[1,1].set_ylim(50000,60000)
# axs[1,1].set_xlim(1.4,2.5)

# axs[2,1].plot(r,Diabatic_Objects["potential"][:,2,3],color="k",alpha=0.25)




# NonAdiabaticCoupling[[2,3]].Rval[2] = 1.90898597152284660

# if Calculation["method"].abinitio_fit == true
#     fit_abinitio()
# elseif Calculation["method"].abinitio_fit == false
#     U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
#     #
#     # fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])

#     # plt.subplots_adjust(wspace=0, hspace=0)

#     # # axs[1,1].set_title("MRCI aug-cc-pVQZ-X2C | occ = 8330, closed = 5110")
#     # for i=1:dim
#     #     if lowercase(Calculation["method"].regularisation) == "potential"
#     #         if i in Calculation["method"].states
#     #             axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
#     #             axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
#     #         end
#     #     else
#     #         for j=i:dim
#     #             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#     #                 axs[1,1].plot(r,Diabatic_Objects[Calculation["method"].regularisation][:,i,j],label="F"*string(i)*string(j))
#     #                 axs[1,1].plot(r,Objects[Calculation["method"].regularisation][:,i,j],"--")
#     #             end
#     #         end
#     #     end
#     # end

#     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], color="green",label= L"$V^{\rm (d)}_1$")
#     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], color="orange",label=L"$V^{\rm (d)}_2$")
#     # # axs[1,1].plot(r,Objects["potential"][:,2,2],color="red","--")
#     # # axs[1,1].plot(r,Objects["potential"][:,3,3],color="blue","--")

#     # # inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25])  # [left, bottom, width, height]

#     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,2,2], color="green")
#     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,3,3], color="orange")
#     # # inset_ax.plot(r,Objects["potential"][:,2,2],color="red","--")
#     # # inset_ax.plot(r,Objects["potential"][:,3,3],color="blue","--")

#     # # inset_ax.set_xlim(1.8, 2.1)  # Zoom into x-range
#     # # inset_ax.set_ylim(50000,63000)  # Zoom into y-range

#     # # axs[1,1].plot([1.8,1.8],[5e4,6.3e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.1],[5e4,6.3e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([1.8,2.1],[5e4,5e4], color="grey"       ,alpha=0.5)
#     # # axs[1,1].plot([1.8,2.1],[6.3e4,6.3e4], color="grey"   ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.668],[5e4,5e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.668],[6.3e4,7.09e4], color="grey",alpha=0.5)

#     # # inset_ax.set_xticklabels([])  # Remove x-axis numbers
#     # # inset_ax.set_yticklabels([])  # Remove y-axis numbers

#     # # # axs[1,1].set_xlabel("Bond Length")
#     # # axs[1,1].set_ylabel(L"Potential Energy, $cm^{-1}$")

#     # for i=1:dim
#     #     for j=i+1:dim
#     #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#     #             # axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
#     #             axs[2,1].plot(r,Objects["nac"][:,i,j]) #,"--",alpha=0.5)
#     #         end
#     #     end
#     # end

#     # # axs[2,1].plot(r,Objects["nac"][:,2,3],"k") #,"--",alpha=0.5)


#     # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
#     # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$")
#     # # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
#     # # axs[1,1].legend()
#     # # axs[2,1].legend()

#     # sf = Calculation["method"].states[end]
#     # Emax = Objects["potential"][end,sf,sf]

#     # si = Calculation["method"].states[1]
#     # Emin = minimum(Objects["potential"][:,si,si])

#     # axs[1,1].set_ylim(Emin,1.1*Emax)
#     # # axs[1,1].set_xlim(1.35,3.5)

#     # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/property_based_pot_guess.png",dpi=300,bbox_inches="tight")
# end


# axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2])
# axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3])

# axs[2,1].plot(r,Diabatic_Objects["potential"][:,2,3])







# Calculation["method"].states = [1,2]


# if Calculation["method"].abinitio_fit == true
#     fit_abinitio()
# elseif Calculation["method"].abinitio_fit == false
#     U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
#     #
#     # fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])

#     # plt.subplots_adjust(wspace=0, hspace=0)

#     # # axs[1,1].set_title("MRCI aug-cc-pVQZ-X2C | occ = 8330, closed = 5110")
#     # for i=1:dim
#     #     if lowercase(Calculation["method"].regularisation) == "potential"
#     #         if i in Calculation["method"].states
#     #             axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
#     #             axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
#     #         end
#     #     else
#     #         for j=i:dim
#     #             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#     #                 axs[1,1].plot(r,Diabatic_Objects[Calculation["method"].regularisation][:,i,j],label="F"*string(i)*string(j))
#     #                 axs[1,1].plot(r,Objects[Calculation["method"].regularisation][:,i,j],"--")
#     #             end
#     #         end
#     #     end
#     # end

#     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], color="green",label= L"$V^{\rm (d)}_1$")
#     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], color="orange",label=L"$V^{\rm (d)}_2$")
#     # # axs[1,1].plot(r,Objects["potential"][:,2,2],color="red","--")
#     # # axs[1,1].plot(r,Objects["potential"][:,3,3],color="blue","--")

#     # # inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25])  # [left, bottom, width, height]

#     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,2,2], color="green")
#     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,3,3], color="orange")
#     # # inset_ax.plot(r,Objects["potential"][:,2,2],color="red","--")
#     # # inset_ax.plot(r,Objects["potential"][:,3,3],color="blue","--")

#     # # inset_ax.set_xlim(1.8, 2.1)  # Zoom into x-range
#     # # inset_ax.set_ylim(50000,63000)  # Zoom into y-range

#     # # axs[1,1].plot([1.8,1.8],[5e4,6.3e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.1],[5e4,6.3e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([1.8,2.1],[5e4,5e4], color="grey"       ,alpha=0.5)
#     # # axs[1,1].plot([1.8,2.1],[6.3e4,6.3e4], color="grey"   ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.668],[5e4,5e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.668],[6.3e4,7.09e4], color="grey",alpha=0.5)

#     # # inset_ax.set_xticklabels([])  # Remove x-axis numbers
#     # # inset_ax.set_yticklabels([])  # Remove y-axis numbers

#     # # # axs[1,1].set_xlabel("Bond Length")
#     # # axs[1,1].set_ylabel(L"Potential Energy, $cm^{-1}$")

#     # for i=1:dim
#     #     for j=i+1:dim
#     #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#     #             # axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
#     #             axs[2,1].plot(r,Objects["nac"][:,i,j]) #,"--",alpha=0.5)
#     #         end
#     #     end
#     # end

#     # # axs[2,1].plot(r,Objects["nac"][:,2,3],"k") #,"--",alpha=0.5)


#     # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
#     # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$")
#     # # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
#     # # axs[1,1].legend()
#     # # axs[2,1].legend()

#     # sf = Calculation["method"].states[end]
#     # Emax = Objects["potential"][end,sf,sf]

#     # si = Calculation["method"].states[1]
#     # Emin = minimum(Objects["potential"][:,si,si])

#     # axs[1,1].set_ylim(Emin,1.1*Emax)
#     # # axs[1,1].set_xlim(1.35,3.5)

#     # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/property_based_pot_guess.png",dpi=300,bbox_inches="tight")
# end
# plt.figure()

# dim=3
# for i=1:dim
#     for j=i:dim
#         plt.plot(r,Diabatic_Objects["potential"][:,i,j])
#         plt.plot(r,Objects["potential"][:,i,j],"--")
#     end
# end

# Calculation["method"].states = [2,3]

# Potential[2].Lval = r
# Potential[2].Rval = Diabatic_Objects["potential"][:,1,1]

# include("Build_Hamiltonian_Matrix.jl")


# if Calculation["method"].abinitio_fit == true
#     fit_abinitio()
# elseif Calculation["method"].abinitio_fit == false
#     U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
#     #
#     # fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])

#     # plt.subplots_adjust(wspace=0, hspace=0)

#     # # axs[1,1].set_title("MRCI aug-cc-pVQZ-X2C | occ = 8330, closed = 5110")
#     # for i=1:dim
#     #     if lowercase(Calculation["method"].regularisation) == "potential"
#     #         if i in Calculation["method"].states
#     #             axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
#     #             axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
#     #         end
#     #     else
#     #         for j=i:dim
#     #             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#     #                 axs[1,1].plot(r,Diabatic_Objects[Calculation["method"].regularisation][:,i,j],label="F"*string(i)*string(j))
#     #                 axs[1,1].plot(r,Objects[Calculation["method"].regularisation][:,i,j],"--")
#     #             end
#     #         end
#     #     end
#     # end

#     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], color="green",label= L"$V^{\rm (d)}_1$")
#     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], color="orange",label=L"$V^{\rm (d)}_2$")
#     # # axs[1,1].plot(r,Objects["potential"][:,2,2],color="red","--")
#     # # axs[1,1].plot(r,Objects["potential"][:,3,3],color="blue","--")

#     # # inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25])  # [left, bottom, width, height]

#     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,2,2], color="green")
#     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,3,3], color="orange")
#     # # inset_ax.plot(r,Objects["potential"][:,2,2],color="red","--")
#     # # inset_ax.plot(r,Objects["potential"][:,3,3],color="blue","--")

#     # # inset_ax.set_xlim(1.8, 2.1)  # Zoom into x-range
#     # # inset_ax.set_ylim(50000,63000)  # Zoom into y-range

#     # # axs[1,1].plot([1.8,1.8],[5e4,6.3e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.1],[5e4,6.3e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([1.8,2.1],[5e4,5e4], color="grey"       ,alpha=0.5)
#     # # axs[1,1].plot([1.8,2.1],[6.3e4,6.3e4], color="grey"   ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.668],[5e4,5e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.668],[6.3e4,7.09e4], color="grey",alpha=0.5)

#     # # inset_ax.set_xticklabels([])  # Remove x-axis numbers
#     # # inset_ax.set_yticklabels([])  # Remove y-axis numbers

#     # # # axs[1,1].set_xlabel("Bond Length")
#     # # axs[1,1].set_ylabel(L"Potential Energy, $cm^{-1}$")

#     # for i=1:dim
#     #     for j=i+1:dim
#     #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#     #             # axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
#     #             axs[2,1].plot(r,Objects["nac"][:,i,j]) #,"--",alpha=0.5)
#     #         end
#     #     end
#     # end

#     # # axs[2,1].plot(r,Objects["nac"][:,2,3],"k") #,"--",alpha=0.5)


#     # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
#     # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$")
#     # # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
#     # # axs[1,1].legend()
#     # # axs[2,1].legend()

#     # sf = Calculation["method"].states[end]
#     # Emax = Objects["potential"][end,sf,sf]

#     # si = Calculation["method"].states[1]
#     # Emin = minimum(Objects["potential"][:,si,si])

#     # axs[1,1].set_ylim(Emin,1.1*Emax)
#     # # axs[1,1].set_xlim(1.35,3.5)

#     # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/property_based_pot_guess.png",dpi=300,bbox_inches="tight")
# end

# plt.figure()

# dim=3
# for i=1:dim
#     for j=i:dim
#         plt.plot(r,Diabatic_Objects["potential"][:,i,j])
#         plt.plot(r,Objects["potential"][:,i,j],"--")
#     end
# end
# 

# dim = 4
# if Calculation["method"].abinitio_fit == true
#     fit_abinitio()
# elseif Calculation["method"].abinitio_fit == false
#     U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
#     # #
#     # fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])

#     # plt.subplots_adjust(wspace=0, hspace=0)

#     # # axs[1,1].set_title("MRCI aug-cc-pVQZ-X2C | occ = 8330, closed = 5110")
#     # for i=1:dim
#     #     if lowercase(Calculation["method"].regularisation) == "potential"
#     #         if i in Calculation["method"].states
#     #             axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
#     #             axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
#     #         end
#     #     else
#     #         for j=i:dim
#     #             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#     #                 axs[1,1].plot(r,Diabatic_Objects[Calculation["method"].regularisation][:,i,j],label="F"*string(i)*string(j))
#     #                 axs[1,1].plot(r,Objects[Calculation["method"].regularisation][:,i,j],"--")
#     #             end
#     #         end
#     #     end
#     # end

#     # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], color="green",label= L"$V^{\rm (d)}_1$")
#     # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], color="orange",label=L"$V^{\rm (d)}_2$")
#     # axs[1,1].plot(r,Objects["potential"][:,2,2],color="red","--")
#     # axs[1,1].plot(r,Objects["potential"][:,3,3],color="blue","--")

#     # inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25])  # [left, bottom, width, height]

#     # inset_ax.plot(r,Diabatic_Objects["potential"][:,2,2], color="green")
#     # inset_ax.plot(r,Diabatic_Objects["potential"][:,3,3], color="orange")
#     # inset_ax.plot(r,Objects["potential"][:,2,2],color="red","--")
#     # inset_ax.plot(r,Objects["potential"][:,3,3],color="blue","--")

#     # inset_ax.set_xlim(1.8, 2.1)  # Zoom into x-range
#     # inset_ax.set_ylim(50000,63000)  # Zoom into y-range

#     # axs[1,1].plot([1.8,1.8],[5e4,6.3e4], color="grey"     ,alpha=0.5)
#     # axs[1,1].plot([2.1,2.1],[5e4,6.3e4], color="grey"     ,alpha=0.5)
#     # axs[1,1].plot([1.8,2.1],[5e4,5e4], color="grey"       ,alpha=0.5)
#     # axs[1,1].plot([1.8,2.1],[6.3e4,6.3e4], color="grey"   ,alpha=0.5)
#     # axs[1,1].plot([2.1,2.668],[5e4,5e4], color="grey"     ,alpha=0.5)
#     # axs[1,1].plot([2.1,2.668],[6.3e4,7.09e4], color="grey",alpha=0.5)

#     # inset_ax.set_xticklabels([])  # Remove x-axis numbers
#     # inset_ax.set_yticklabels([])  # Remove y-axis numbers

#     # # axs[1,1].set_xlabel("Bond Length")
#     # # axs[1,1].set_ylabel(L"Potential Energy, $cm^{-1}$")

#     # for i=1:dim
#     #     for j=i+1:dim
#     #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#     #             # axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
#     #             axs[2,1].plot(r,Objects["nac"][:,i,j]) #,"--",alpha=0.5)
#     #         end
#     #     end
#     # end

#     # # axs[2,1].plot(r,Objects["nac"][:,2,3],"k") #,"--",alpha=0.5)


#     # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
#     # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$")
#     # # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
#     # # axs[1,1].legend()
#     # # axs[2,1].legend()

#     # sf = Calculation["method"].states[end]
#     # Emax = Objects["potential"][end,sf,sf]

#     # si = Calculation["method"].states[1]
#     # Emin = minimum(Objects["potential"][:,si,si])

#     # axs[1,1].set_ylim(Emin,1.1*Emax)
#     # axs[1,1].set_xlim(1.35,3.5)

#     # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/property_based_pot_guess.png",dpi=300,bbox_inches="tight")
# end

# Calculation["method"].diabatisation = "forward-evolution"

# U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))

# fig, axs = plt.subplots(2, 2,sharex=true)
# #
# axs[1,1].tick_params(axis="both", labelsize=8)
# axs[1,2].tick_params(axis="both", labelsize=8)
# axs[2,1].tick_params(axis="both", labelsize=8)
# axs[2,2].tick_params(axis="both", labelsize=8)
# lw = 0.9
# #
# axs[1,1].set_title("Forward Evolution",fontsize=8)

# axs[1,1].plot(r,Diabatic_Objects["potential"][:,1,1], linewidth=lw,color="black",label = L"V$^{\rm(d)}_{1}$")
# axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], linewidth=lw,color="grey",label = L"V$^{\rm(d)}_{2}$")
# axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], linewidth=lw,color="orange",label = L"V$^{\rm(d)}_{3}$")
# axs[1,1].plot(r,Diabatic_Objects["potential"][:,4,4], linewidth=lw,color="brown",label = L"V$^{\rm(d)}_{4}$")

# axs[1,1].plot(r,Objects["potential"][:,1,1], "--", alpha=0.5,linewidth=lw, color="red",    label = L"$C^{2}\Sigma^{+}+\epsilon K_{11}$")
# axs[1,1].plot(r,Objects["potential"][:,2,2], "--", alpha=0.5,linewidth=lw, color="green",   label = L"$(2)^{2}\Sigma^{+}+\epsilon K_{22}$")
# axs[1,1].plot(r,Objects["potential"][:,3,3], "--", alpha=0.5,linewidth=lw, color="blue",   label = L"$(3)^{2}\Sigma^{+}+\epsilon K_{33}$")
# axs[1,1].plot(r,Objects["potential"][:,4,4], "--", alpha=0.5,linewidth=lw, color="magenta",label = L"$(4)^{2}\Sigma^{+}+\epsilon K_{44}$")

# axs[2,1].plot(r,Diabatic_Objects["potential"][:,2,3],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{23}$")
# axs[2,1].plot(r,Diabatic_Objects["potential"][:,1,3],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{13}$")
# axs[2,1].plot(r,Diabatic_Objects["potential"][:,1,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{14}$")
# axs[2,1].plot(r,Diabatic_Objects["potential"][:,2,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{24}$")
# axs[2,1].plot(r,Diabatic_Objects["potential"][:,3,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{34}$")
# axs[2,1].plot(r,Diabatic_Objects["potential"][:,1,2],linewidth=lw,color="red",label=L"W$^{\rm(1)}_{12}$")



# Calculation["method"].diabatisation = "backward-evolution"

# U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))


# axs[1,2].set_title("Backward Evolution",fontsize=8)

# axs[1,2].plot(r,Diabatic_Objects["potential"][:,1,1], linewidth=lw,color="black",label = L"V$^{\rm(d)}_{1}$")
# axs[1,2].plot(r,Diabatic_Objects["potential"][:,2,2], linewidth=lw,color="grey",label = L"V$^{\rm(d)}_{2}$")
# axs[1,2].plot(r,Diabatic_Objects["potential"][:,3,3], linewidth=lw,color="orange",label = L"V$^{\rm(d)}_{3}$")
# axs[1,2].plot(r,Diabatic_Objects["potential"][:,4,4], linewidth=lw,color="brown",label = L"V$^{\rm(d)}_{4}$")

# axs[1,2].plot(r,Objects["potential"][:,1,1], "--", linewidth=lw, alpha=0.5,color="red",    label = L"$C^{2}\Sigma^{+}+\epsilon K_{11}$")
# axs[1,2].plot(r,Objects["potential"][:,2,2], "--", linewidth=lw, alpha=0.5,color="green",   label = L"$(2)^{2}\Sigma^{+}+\epsilon K_{22}$")
# axs[1,2].plot(r,Objects["potential"][:,3,3], "--", linewidth=lw, alpha=0.5,color="blue",   label = L"$(3)^{2}\Sigma^{+}+\epsilon K_{33}$")
# axs[1,2].plot(r,Objects["potential"][:,4,4], "--", linewidth=lw, alpha=0.5,color="magenta",label = L"$(4)^{2}\Sigma^{+}+\epsilon K_{44}$")

# axs[2,2].plot(r,Diabatic_Objects["potential"][:,2,3],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{23}$")
# axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,3],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{13}$")
# axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{14}$")
# axs[2,2].plot(r,Diabatic_Objects["potential"][:,2,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{24}$")
# axs[2,2].plot(r,Diabatic_Objects["potential"][:,3,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{34}$")
# axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,2],linewidth=lw,color="red",label=L"W$^{\rm(1)}_{12}$")


# #
# axs[2,1].set_xlabel(L"Bond Length R, $\mathrm{\AA}$",fontsize=8)
# axs[2,2].set_xlabel(L"Bond Length R, $\mathrm{\AA}$",fontsize=8)
# #
# axs[2,1].set_ylabel(L"DC, $\rm cm^{-1}$",fontsize=8)
# #
# axs[1,1].set_ylabel(L"Potential, cm$^{-1}$",fontsize=8)
# #
# axs[1,1].set_ylim(25000,120000)
# axs[1,2].set_ylim(25000,120000)
# axs[2,1].set_ylim(-10100,21000)
# axs[2,2].set_ylim(-10100,21000)
# axs[1,1].set_xlim(r[1],6)
# axs[2,1].set_xlim(r[1],6)
# axs[1,2].set_xlim(r[1],6)
# axs[2,2].set_xlim(r[1],6)

# axs[1,2].set_yticklabels([])
# axs[2,2].set_yticklabels([])

# plt.subplots_adjust(hspace=0.0,wspace=0.1)
# fig.set_size_inches(6.69, 3.1)

# plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/N-level/Evolution/CH/plots/paper/CH_forward_backward_evo_diff_PECs_DCs_simplified.png",bbox_inches="tight",dpi=600)









# plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/CH_scale_NAC.png",bbox_inches="tight",dpi=600)
# plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/YO_shift_NAC_2.png",bbox_inches="tight",dpi=600)

# plt.figure()

# for i=1:dim
#     for j=i:dim
#         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#             plt.plot(r,Objects["nac"][:,i,j])
#             # plt.plot(r,[Objects["regularised_nac"][idx][i,j] for idx=1:lastindex(r)],"--")
#         end
#     end
# end

# save_diabatisation(Objects, Diabatic_Objects, lowercase(Calculation["method"].diabatisation), input_properties, fname, special_name = "4paper")



# fig, axs = plt.subplots(2, 2,sharex=true)
# axs[1,1].tick_params(axis="both", labelsize=8)
# axs[1,2].tick_params(axis="both", labelsize=8)
# axs[2,1].tick_params(axis="both", labelsize=8)
# axs[2,2].tick_params(axis="both", labelsize=8)
# lw = 0.9
# #
# ## diabats
# axs[1,2].text(5.2,110000, L"\rm CH",fontsize=10,color="black")

# axs[1,2].plot(r,Diabatic_Objects["potential"][:,1,1], linewidth=lw, color="black",label = L"V$^{\rm(d)}_{1}$")
# axs[1,2].plot(r,Diabatic_Objects["potential"][:,2,2], linewidth=lw, color="grey",label = L"V$^{\rm(d)}_{2}$")
# axs[1,2].plot(r,Diabatic_Objects["potential"][:,3,3], linewidth=lw, color="orange",label = L"V$^{\rm(d)}_{3}$")
# axs[1,2].plot(r,Diabatic_Objects["potential"][:,4,4], linewidth=lw, color="blue",label = L"V$^{\rm(d)}_{4}$")


# axs[1,2].text(1.83,59100,  L"V$^{\rm(d)}_{3}$",fontsize=8,color="orange")
# axs[1,2].text(1.58 ,31000,  L"V$^{\rm(d)}_{1}$",fontsize=8,color="black")
# axs[1,2].text(3.15,108000,  L"V$^{\rm(d)}_{2}$",fontsize=8,color="grey")
# axs[1,2].text(2.87,86300,  L"V$^{\rm(d)}_{4}$",fontsize=8,color="blue")

# axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,2], linewidth=lw, color="blue", label=L"$\mathcal{D}_{12}$")
# axs[2,2].plot(r,Diabatic_Objects["potential"][:,2,3], linewidth=lw, color="red",  label=L"$\mathcal{D}_{23}$")
# axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,3], linewidth=lw, color="green",label=L"$\mathcal{D}_{13}$")
# axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,4], linewidth=lw, color="black",label=L"$\mathcal{D}_{14}$")
# axs[2,2].plot(r,Diabatic_Objects["potential"][:,2,4], linewidth=lw, color="orange",label=L"$\mathcal{D}_{24}$")
# axs[2,2].plot(r,Diabatic_Objects["potential"][:,3,4], linewidth=lw, color="purple",label=L"$\mathcal{D}_{34}$")


# axs[2,2].text(1.25,3000,  L"$V^{\rm(d)}_{12}$",fontsize=8,color="blue")
# axs[2,2].text(2.33, 10000, L"$V^{\rm(d)}_{23}$",fontsize=8,color="red")
# axs[2,2].text(3,-4500 , L"$V^{\rm(d)}_{13}$",fontsize=8,color="green")
# axs[2,2].text(4.5,3000, L"$V^{\rm(d)}_{14}$",fontsize=8,color="black")
# axs[2,2].text(3.2 ,4300,  L"$V^{\rm(d)}_{24}$",fontsize=8,color="orange")
# axs[2,2].text(2.25,-2500,  L"$V^{\rm(d)}_{34}$",fontsize=8,color="purple")

# #
# ## adiabats

# atom1, atom2 = Calculation["method"].atoms
# kinetic_factor =  KE_factor(atom1, atom2)


# axs[1,1].plot(r,Objects["potential"][:,1,1] .+ (kinetic_factor .* [k[1,1] for k in Objects["K_matrix"]]), linewidth=lw, color="red",    label = L"$C^{2}\Sigma^{+}+\epsilon K_{11}$")
# axs[1,1].plot(r,Objects["potential"][:,2,2] .+ (kinetic_factor .* [k[2,2] for k in Objects["K_matrix"]]), linewidth=lw, color="green",   label = L"$(2)^{2}\Sigma^{+}+\epsilon K_{22}$")
# axs[1,1].plot(r,Objects["potential"][:,3,3] .+ (kinetic_factor .* [k[3,3] for k in Objects["K_matrix"]]), linewidth=lw, color="blue",   label = L"$(3)^{2}\Sigma^{+}+\epsilon K_{33}$")
# axs[1,1].plot(r,Objects["potential"][:,4,4] .+ (kinetic_factor .* [k[4,4] for k in Objects["K_matrix"]]), linewidth=lw, color="magenta",label = L"$(4)^{2}\Sigma^{+}+\epsilon K_{44}$")


# axs[1,1].text(3.5,33000,  L"$V_{C^{2}\Sigma^{+}}+\epsilon K_{11}$",fontsize=8,color="red")
# axs[1,1].text(3.5,44000,  L"$V_{2^{2}\Sigma^{+}}+\epsilon K_{22}$",fontsize=8,color="green")
# axs[1,1].text(3.5,72000,  L"$V_{3^{2}\Sigma^{+}}+\epsilon K_{33}$",fontsize=8,color="blue")
# axs[1,1].text(3.5,110000, L"$V_{4^{2}\Sigma^{+}}+\epsilon K_{44}$",fontsize=8,color="magenta")


# axs[2,1].plot(r,[w[1,2] for w in Objects["regularised_nac"]], linewidth=lw,color="blue",label=L"W$^{\rm(1)}_{12}$")
# axs[2,1].plot(r,[w[2,3] for w in Objects["regularised_nac"]], linewidth=lw,color="red",label=L"W$^{\rm(1)}_{23}$")
# axs[2,1].plot(r,[w[1,3] for w in Objects["regularised_nac"]], linewidth=lw,color="green",label=L"W$^{\rm(1)}_{13}$")
# axs[2,1].plot(r,[w[1,4] for w in Objects["regularised_nac"]], linewidth=lw,color="black",label=L"W$^{\rm(1)}_{14}$")
# axs[2,1].plot(r,[w[2,4] for w in Objects["regularised_nac"]], linewidth=lw,color="orange",label=L"W$^{\rm(1)}_{24}$")
# axs[2,1].plot(r,[w[3,4] for w in Objects["regularised_nac"]], linewidth=lw,color="purple",label=L"W$^{\rm(1)}_{34}$")

# axs[2,1].text(2.95, 1.4  , L"$\langle C^{2}\Sigma^{+}|\frac{d}{dr}|2^{2}\Sigma^{+}\rangle$",fontsize=8,color="blue")
# axs[2,1].text(2.95, 1.05  ,L"$\langle 2^{2}\Sigma^{+}|\frac{d}{dr}|3^{2}\Sigma^{+}\rangle$",fontsize=8,color="red")
# axs[2,1].text(2.95, 0.35 , L"$\langle C^{2}\Sigma^{+}|\frac{d}{dr}|3^{2}\Sigma^{+}\rangle$",fontsize=8,color="green")
# axs[2,1].text(2.95, -0.25 , L"$\langle C^{2}\Sigma^{+}|\frac{d}{dr}|4^{2}\Sigma^{+}\rangle$",fontsize=8,color="black")
# axs[2,1].text(2.95, 0.7  , L"$\langle 2^{2}\Sigma^{+}|\frac{d}{dr}|4^{2}\Sigma^{+}\rangle$",fontsize=8,color="orange")
# axs[2,1].text(2.95, 1.75 , L"$\langle 3^{2}\Sigma^{+}|\frac{d}{dr}|4^{2}\Sigma^{+}\rangle$",fontsize=8,color="purple")
# #
# axs[2,1].set_xlabel(L"Bond Length R, $\mathrm{\AA}$",fontsize=8)
# axs[2,2].set_xlabel(L"Bond Length R, $\mathrm{\AA}$",fontsize=8)
# #
# axs[2,1].set_ylabel("Coupling",fontsize=8)
# #
# axs[1,1].set_ylabel(L"Potential Energy, cm$^{-1}$",fontsize=8)
# #
# axs[1,1].set_ylim(25000,120000)
# axs[1,2].set_ylim(25000,120000)
# axs[1,1].set_xlim(r[1],6)
# axs[2,1].set_xlim(r[1],6)

# plt.subplots_adjust(hspace=0.0,wspace=0.4)
# fig.set_size_inches(4.5,5)

# axs[1,1].tick_params(axis="both", labelsize=8)
# axs[2,1].tick_params(axis="both", labelsize=8)
# axs[1,2].tick_params(axis="both", labelsize=8)
# axs[2,2].tick_params(axis="both", labelsize=8)

# plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/wavefunction_comparer/CH/plots/CH_diabatisation_adi-dia_4x5.png",bbox_inches="tight",dpi=600)





#
######## CO

# Va2 = deepcopy(Potential[2].Rval)
# rVa2 = deepcopy(Potential[2].Lval)


# QD = deepcopy(Diabatic_Objects["potential"][:,3,3])
# R = deepcopy(r)

# Calculation["method"].states = [1,2]

# Potential[2].Lval = r
# Potential[2].Rval = QD

# if Calculation["method"].abinitio_fit == true
#     fit_abinitio()
# elseif Calculation["method"].abinitio_fit == false
#     U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
#     #
#     # fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])

#     # plt.subplots_adjust(wspace=0, hspace=0)

#     # # axs[1,1].set_title("MRCI aug-cc-pVQZ-X2C | occ = 8330, closed = 5110")
#     # for i=1:dim
#     #     if lowercase(Calculation["method"].regularisation) == "potential"
#     #         if i in Calculation["method"].states
#     #             axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
#     #             axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
#     #         end
#     #     else
#     #         for j=i:dim
#     #             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#     #                 axs[1,1].plot(r,Diabatic_Objects[Calculation["method"].regularisation][:,i,j],label="F"*string(i)*string(j))
#     #                 axs[1,1].plot(r,Objects[Calculation["method"].regularisation][:,i,j],"--")
#     #             end
#     #         end
#     #     end
#     # end

#     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], color="green",label= L"$V^{\rm (d)}_1$")
#     # # axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], color="orange",label=L"$V^{\rm (d)}_2$")
#     # # axs[1,1].plot(r,Objects["potential"][:,2,2],color="red","--")
#     # # axs[1,1].plot(r,Objects["potential"][:,3,3],color="blue","--")

#     # # inset_ax = fig.add_axes([0.6, 0.6, 0.25, 0.25])  # [left, bottom, width, height]

#     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,2,2], color="green")
#     # # inset_ax.plot(r,Diabatic_Objects["potential"][:,3,3], color="orange")
#     # # inset_ax.plot(r,Objects["potential"][:,2,2],color="red","--")
#     # # inset_ax.plot(r,Objects["potential"][:,3,3],color="blue","--")

#     # # inset_ax.set_xlim(1.8, 2.1)  # Zoom into x-range
#     # # inset_ax.set_ylim(50000,63000)  # Zoom into y-range

#     # # axs[1,1].plot([1.8,1.8],[5e4,6.3e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.1],[5e4,6.3e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([1.8,2.1],[5e4,5e4], color="grey"       ,alpha=0.5)
#     # # axs[1,1].plot([1.8,2.1],[6.3e4,6.3e4], color="grey"   ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.668],[5e4,5e4], color="grey"     ,alpha=0.5)
#     # # axs[1,1].plot([2.1,2.668],[6.3e4,7.09e4], color="grey",alpha=0.5)

#     # # inset_ax.set_xticklabels([])  # Remove x-axis numbers
#     # # inset_ax.set_yticklabels([])  # Remove y-axis numbers

#     # # # axs[1,1].set_xlabel("Bond Length")
#     # # axs[1,1].set_ylabel(L"Potential Energy, $cm^{-1}$")

#     # for i=1:dim
#     #     for j=i+1:dim
#     #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#     #             # axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
#     #             axs[2,1].plot(r,Objects["nac"][:,i,j]) #,"--",alpha=0.5)
#     #         end
#     #     end
#     # end

#     # # axs[2,1].plot(r,Objects["nac"][:,2,3],"k") #,"--",alpha=0.5)


#     # axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
#     # axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$")
#     # # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
#     # # axs[1,1].legend()
#     # # axs[2,1].legend()

#     # sf = Calculation["method"].states[end]
#     # Emax = Objects["potential"][end,sf,sf]

#     # si = Calculation["method"].states[1]
#     # Emin = minimum(Objects["potential"][:,si,si])

#     # axs[1,1].set_ylim(Emin,1.1*Emax)
#     # # axs[1,1].set_xlim(1.35,3.5)

#     # plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/Thesis/AtDT_method_comp/property_based_pot_guess.png",dpi=300,bbox_inches="tight")
# end

# rep = deepcopy(Diabatic_Objects["potential"][:,2,2])



# f = NonAdiabaticCoupling[[1,2]].type
# p =  NonAdiabaticCoupling[[1,2]].Rval
# #
# ##
# println("CATS")
# mixing_angle_ij = compute_mixing_angle(r, f, [1,2], p)
# #
# c = cos.(mixing_angle_ij)
# s = sin.(mixing_angle_ij)
# #
# U = zeros(Float64,length(r),2,2)
# U[:,1,1] .= c
# U[:,2,2] .= c
# U[:,1,2] .= -s
# U[:,2,1] .= s

# Va = zeros(Float64,length(r),2,2)
# Va[:,1,1] .= Objects["potential"][:,1,1]
# Va[:,2,2] .= QD

# Vd = map(idx -> U[idx,:,:]' * Va[idx,:,:] * U[idx,:,:], collect(1:lastindex(r)))

# rep = [v[2,2] for v in Vd]


# Calculation["grid"].npoints = 3001
# include("Build_Hamiltonian_Matrix.jl")

# fig, axs = plt.subplots(2,2,figsize=[6,4],sharex=true)
# plt.subplots_adjust(hspace=0, wspace=0.65)  # Remove vertical space


# # plt.plot(Potential[5].Lval,Potential[5].Rval)
# # plt.plot(Potential[6].Lval,Potential[6].Rval)
# # plt.plot(Potential[7].Lval,Potential[7].Rval)

# axs[1,1].plot(Potential[1].Lval,Potential[1].Rval,"blue")
# axs[1,1].plot(Potential[3].Lval,Potential[3].Rval,"orange")
# axs[1,1].plot(rVa2,Va2,"limegreen")
# # axs[1,1].plot(Potential[4].Lval,Potential[4].Rval)
# axs[1,1].plot(R,QD,color="red","--")
# axs[1,1].set_ylabel(L"Potential Energy, cm$^{-1}$")
# axs[1,1].text(1.4,92000,L"$B^1\Sigma^+$",fontsize=9,color="blue")
# axs[1,1].text(1.06,97000,L"$C^1\Sigma^+$",fontsize=9,color="limegreen")
# axs[1,1].text(1.1,113000,L"$(IV)^1\Sigma^+$",fontsize=9,color="orange")
# axs[1,1].text(1.4,100000,L"$QD^1\Sigma^+$",fontsize=9,color="red")


# axs[1,1].set_ylim(85000,120000)


# global r
# r = LinRange(Calculation["grid"].range[1],
#              Calculation["grid"].range[2],
#              Calculation["grid"].npoints)  

# axs[2,1].plot(r,Objects["nac"][:,1,2],color="k")
# axs[2,1].plot(r,Objects["nac"][:,2,3],color="red")
# axs[2,1].set_ylabel(L"NAC, $\rm \AA^{-1}$" )
# axs[2,1].set_xlabel(L"Bond Length, $\rm \AA$")
# axs[2,1].text(0.95,45,L"$\langle C^1\Sigma^+| d/dr | (IV)^1\Sigma^+ \rangle$",color="red",fontsize=9)
# axs[2,1].text(1.1,40,L"$\langle B^1\Sigma^+| d/dr | QD^1\Sigma^+ \rangle$",color="k",fontsize=9)
# axs[2,1].set_ylim(0,50)

# axs[1,2].plot(R,rep,"orange",alpha=0.25)
# axs[1,2].plot(Potential[5].Lval,Potential[5].Rval,"blue")
# axs[1,2].plot(Potential[6].Lval,Potential[6].Rval,"green")
# axs[1,2].set_ylabel(L"Potential Energy, cm$^{-1}$")
# axs[1,2].text(1.175,105000,L"$C^1\Sigma^{+\rm(d)}$",fontsize=9,color="green")
# axs[1,2].text(1.4,95000,L"$B^1\Sigma^{+\rm(d)}$",fontsize=9,color="blue")
# axs[1,2].text(1.1,115000,L"$(3)^1\Sigma^{+\rm(d)}$",fontsize=9,color="orange")



# axs[2,2].plot(Potential[8].Lval,Potential[8].Rval.*219474.6313708000,color="black")
# axs[2,2].plot(Potential[9].Lval,Potential[9].Rval,color="red")
# axs[2,2].set_ylabel(L"DC, cm$^{-1}$")
# axs[2,2].set_xlabel(L"Bond Length, $\rm \AA$")
# axs[2,2].set_ylim(0,4000)

# axs[2,2].text(0.95,3550,L"$\langle C^1\Sigma^{+\rm(d)}| V^{\rm(d)} | QD^1\Sigma^+ \rangle$",color="red",fontsize=9)
# axs[2,2].text(0.95,3000,L"$\langle B^1\Sigma^{+\rm(d)}| V^{\rm(d)} | (3)^1\Sigma^{+ \rm (d)} \rangle$",color="k",fontsize=9)
# # axs[2,2].text(1.4,500,L"$\times 50000$",color="k",fontsize=8)


# axs[1,2].set_xlim(0.9,1.6)
# axs[1,2].set_ylim(85000,120000)

# plt.savefig("./CO_2state_diabatisation_2.png", bbox_inches="tight",dpi=300)







