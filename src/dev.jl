
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
fname = "/Users/ryanbrady/Documents/PhD/Work/DIABATISATION/DIABATOM_PRO_PACKAGE/github_dev/DIABATOMPRO/Supplementary/KH/KH.inp"
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
#     U, dU, UdU, K_Matrix, diabatic_basis, Diabatic_Objects, input_properties = run_diabatiser(lowercase(Calculation["method"].diabatisation))
#     #
#     fig, axs = plt.subplots(2,1,sharex=true,figsize=[8,8])

#     plt.subplots_adjust(wspace=0, hspace=0)

#     # axs[1,1].set_title("MRCI aug-cc-pVQZ-X2C | occ = 8330, closed = 5110")
#     for i=1:dim
#         if i in Calculation["method"].states
#             axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
#             axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
#         end
#     end

#     # axs[1,1].set_xlabel("Bond Length")
#     axs[1,1].set_ylabel("Potential, cm-1")

#     for i=1:dim
#         for j=i+1:dim
#             if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#                 axs[2,1].plot(r,Objects["regularised_nac"][:,i,j],label="<"*string(i)*"| d/dr |"*string(j)*">")
#                 axs[2,1].plot(r,Objects["nac"][:,i,j],"--",alpha=0.5)
#             end
#         end
#     end

#     axs[2,1].set_xlabel("Bond Length")
#     axs[2,1].set_ylabel("NAC, 1/Ang")
#     # plt.plot(r,Diabatic_Objects["potential"][:,3,3])
#     axs[1,1].legend()
#     axs[2,1].legend()

#     sf = Calculation["method"].states[end]
#     Emax = Objects["potential"][end,sf,sf]

#     si = Calculation["method"].states[1]
#     Emin = minimum(Objects["potential"][:,si,si])

#     axs[1,1].set_ylim(Emin,1.1*Emax)
# end
#
## make a new Hamiltonian with different representations
#
## save diabatisation?
