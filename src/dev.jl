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
fname =  "/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/DIABATOM_PRO_PACKAGE/github_dev/DIABATOMPRO/Supplementary/CH/DIABATOM-PRO_CH.inp" #"/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/DIABATOM_PRO_PACKAGE/github_dev/DIABATOMPRO/Supplementary/CH/2Pi/CH_doublet_pi.inp" "/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/DIABATOM_PRO_PACKAGE/github_dev/DIABATOMPRO/Supplementary/SO/SO.inp"
read_file(fname)
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
#
## initialise diabatom object
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
#
## run the diabatiser
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

# for i=1:dim
#     for j=i:dim
#         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#             plt.plot(r,Diabatic_Objects["dipole"][:,i,j])
#             plt.plot(r,Objects["dipole"][:,i,j],"--")
#         end
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

Calculation["method"].diabatisation = "forward-evolution"

U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))

fig, axs = plt.subplots(2, 2,sharex=true)
#
axs[1,1].tick_params(axis="both", labelsize=8)
axs[1,2].tick_params(axis="both", labelsize=8)
axs[2,1].tick_params(axis="both", labelsize=8)
axs[2,2].tick_params(axis="both", labelsize=8)
lw = 0.9
#
axs[1,1].set_title("Forward Evolution",fontsize=8)

axs[1,1].plot(r,Diabatic_Objects["potential"][:,1,1], linewidth=lw,color="black",label = L"V$^{\rm(d)}_{1}$")
axs[1,1].plot(r,Diabatic_Objects["potential"][:,2,2], linewidth=lw,color="grey",label = L"V$^{\rm(d)}_{2}$")
axs[1,1].plot(r,Diabatic_Objects["potential"][:,3,3], linewidth=lw,color="orange",label = L"V$^{\rm(d)}_{3}$")
axs[1,1].plot(r,Diabatic_Objects["potential"][:,4,4], linewidth=lw,color="brown",label = L"V$^{\rm(d)}_{4}$")

axs[1,1].plot(r,Objects["potential"][:,1,1], "--", linewidth=lw, color="red",    label = L"$C^{2}\Sigma^{+}+\epsilon K_{11}$")
axs[1,1].plot(r,Objects["potential"][:,2,2], "--", linewidth=lw, color="green",   label = L"$(2)^{2}\Sigma^{+}+\epsilon K_{22}$")
axs[1,1].plot(r,Objects["potential"][:,3,3], "--", linewidth=lw, color="blue",   label = L"$(3)^{2}\Sigma^{+}+\epsilon K_{33}$")
axs[1,1].plot(r,Objects["potential"][:,4,4], "--", linewidth=lw, color="magenta",label = L"$(4)^{2}\Sigma^{+}+\epsilon K_{44}$")

axs[2,1].plot(r,Diabatic_Objects["potential"][:,2,3],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{23}$")
axs[2,1].plot(r,Diabatic_Objects["potential"][:,1,3],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{13}$")
axs[2,1].plot(r,Diabatic_Objects["potential"][:,1,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{14}$")
axs[2,1].plot(r,Diabatic_Objects["potential"][:,2,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{24}$")
axs[2,1].plot(r,Diabatic_Objects["potential"][:,3,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{34}$")
axs[2,1].plot(r,Diabatic_Objects["potential"][:,1,2],linewidth=lw,color="red",label=L"W$^{\rm(1)}_{12}$")



Calculation["method"].diabatisation = "backward-evolution"

U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))


axs[1,2].set_title("Backward Evolution",fontsize=8)

axs[1,2].plot(r,Diabatic_Objects["potential"][:,1,1], linewidth=lw,color="black",label = L"V$^{\rm(d)}_{1}$")
axs[1,2].plot(r,Diabatic_Objects["potential"][:,2,2], linewidth=lw,color="grey",label = L"V$^{\rm(d)}_{2}$")
axs[1,2].plot(r,Diabatic_Objects["potential"][:,3,3], linewidth=lw,color="orange",label = L"V$^{\rm(d)}_{3}$")
axs[1,2].plot(r,Diabatic_Objects["potential"][:,4,4], linewidth=lw,color="brown",label = L"V$^{\rm(d)}_{4}$")

axs[1,2].plot(r,Objects["potential"][:,1,1], "--", linewidth=lw, color="red",    label = L"$C^{2}\Sigma^{+}+\epsilon K_{11}$")
axs[1,2].plot(r,Objects["potential"][:,2,2], "--", linewidth=lw, color="green",   label = L"$(2)^{2}\Sigma^{+}+\epsilon K_{22}$")
axs[1,2].plot(r,Objects["potential"][:,3,3], "--", linewidth=lw, color="blue",   label = L"$(3)^{2}\Sigma^{+}+\epsilon K_{33}$")
axs[1,2].plot(r,Objects["potential"][:,4,4], "--", linewidth=lw, color="magenta",label = L"$(4)^{2}\Sigma^{+}+\epsilon K_{44}$")

axs[2,2].plot(r,Diabatic_Objects["potential"][:,2,3],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{23}$")
axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,3],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{13}$")
axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{14}$")
axs[2,2].plot(r,Diabatic_Objects["potential"][:,2,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{24}$")
axs[2,2].plot(r,Diabatic_Objects["potential"][:,3,4],linewidth=lw,color="black",alpha=0.25,label=L"W$^{\rm(1)}_{34}$")
axs[2,2].plot(r,Diabatic_Objects["potential"][:,1,2],linewidth=lw,color="red",label=L"W$^{\rm(1)}_{12}$")


#
axs[2,1].set_xlabel(L"Bond Length R, $\mathrm{\AA}$",fontsize=8)
axs[2,2].set_xlabel(L"Bond Length R, $\mathrm{\AA}$",fontsize=8)
#
axs[2,1].set_ylabel(L"DC, $\rm cm^{-1}$",fontsize=8)
#
axs[1,1].set_ylabel(L"Potential, cm$^{-1}$",fontsize=8)
#
axs[1,1].set_ylim(25000,120000)
axs[1,2].set_ylim(25000,120000)
axs[2,1].set_ylim(-10100,21000)
axs[2,2].set_ylim(-10100,21000)
axs[1,1].set_xlim(r[1],6)
axs[2,1].set_xlim(r[1],6)
axs[1,2].set_xlim(r[1],6)
axs[2,2].set_xlim(r[1],6)

axs[1,2].set_yticklabels([])
axs[2,2].set_yticklabels([])

plt.subplots_adjust(hspace=0.0,wspace=0.1)
fig.set_size_inches(6.69, 3.1)

plt.savefig("/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/N-level/Evolution/CH/plots/paper/CH_forward_backward_evo_diff_PECs_DCs_simplified.png",bbox_inches="tight",dpi=600)









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
