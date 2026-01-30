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
read_file("../Supplementary/CN_4_states.inp")
#~~~~~~~~~~~~~~~~~~~~~~~~~ RUN HAMILTONIAN BUILDER ~~~~~~~~~~~~~~~~~~~~~~~~#
include("Build_Hamiltonian_Matrix.jl")
Objects = build_hamiltonian_objects(Calculation, 
                                               Potential, 
                                               SpinOrbit, 
                                               EAMC, 
                                               Dipole, 
                                               NonAdiabaticCoupling)
#~~~~~~~~~~~~~~~~~~~~~~~ RUN DIABATISATION PIPELINES ~~~~~~~~~~~~~~~~~~~~~~#
#
## define the bond length vector
dim = length(Potential)
global r
r = LinRange(Calculation["grid"].range[1],
             Calculation["grid"].range[2],
             Calculation["grid"].npoints)  

r  = collect(r)


# fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])
# plt.subplots_adjust(wspace=0, hspace=0)

# states = [1,2,3,6]
# plt.figure()
# for i in states
#     for j in states
#         if i <= j
#             if [i,j] in keys(SpinOrbit)
#                 #

#                 SO = SpinOrbit[[i,j]] #Hamiltonian[("spin-orbit",[i,j])]
#                 #
#                 Li = Vi.lambda
#                 #
#                 Si = Vi.mult
#                 #
#                 Lj = Vj.lambda
#                 #
#                 Sj = Vj.mult
#                 #
#                 if (Si == Sj) & (Li == Lj)
#                     axs[1,1].plot(r,Objects["spin-orbit"][:,i,i],label="<" * string(i) * "|SOZ|" * string(j)*">")
#                 else
#                     axs[1,1].plot(r,Objects["spin-orbit"][:,i,i],label="<" * string(i) * "|SOX|" * string(j)*">")
#                 end
#             end
#         end
#     end
# end
# plt.legend()


# compute vibronic energies and wavefunctions for a non-rotating molecule if
# the 'vibronic_solver' key is in Calculation
# run the diabatiser
if Calculation["method"].abinitio_fit == true
    fit_abinitio()
elseif Calculation["method"].abinitio_fit == false
    U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
   

    # for key in keys(Hamiltonian)
    #     obj = Hamiltonian[key]
    #     #
    #     obj.Rval =
    

    # Calculation["method"].states = 
    
    # fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])
    # plt.subplots_adjust(wspace=0, hspace=0)
    # # for i=1:dim
    # #     if i in Calculation["method"].states
    # #         axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
    # #         axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
    # #     end
    # # end
    # #
    # for j=16:17
    #     axs[1,1].plot(r,Diabatic_Objects["lx"][:,13,j],label="dia 13|"*string(j))
    #     axs[1,1].plot(r,Objects["lx"][:,13,j],"--",label="adi 13|"*string(j))
    #     #
    #     axs[1,1].plot(r,Diabatic_Objects["lx"][:,14,j],label="dia 14|"*string(j))
    #     axs[1,1].plot(r,Objects["lx"][:,14,j],"--",label="adi 14|"*string(j))
    # end
    # axs[1,1].legend()
    # #
    # for i=1:dim
    #     for j=i+1:dim
    #         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
    #             axs[2,1].plot(r,Objects["nac"][:,i,j])
    #             # if (lowercase(Calculation["method"].diabatisation) == "evolution")
    #             #     axs[2,1].plot(r,[Objects["regularised_nac"][idx][i,j] for idx=1:lastindex(r)],"--",label="reg $i $j")
    #             # end
    #         end
    #     end
    # end
    # axs[2,1].set_xlabel("Bond Length, Angstrom")
    # axs[2,1].set_ylabel("NAC, 1/Angstrom")
    # axs[1,1].set_ylabel("Potential, cm-1")
    # axs[2,1].legend()

    # axs[1,1].set_ylim(50000,85000)
    # axs[1,1].set_xlim(1,1.5)


    # plt.figure()
    # plt.plot(r,map(x -> Objects["K_matrix"][x][1,1], collect(1:lastindex(r))))
    # plt.xlabel("Bond Length, Angstroms")
    # plt.ylabel("DBOC, cm-1")

    # plt.figure()
    # plt.plot(r,Diabatic_Objects["potential"][:,1,2])
    # plt.xlabel("Bond Length, Angstroms")
    # plt.ylabel("Diabatic (potential) coupling, cm-1")

end


fig, axs = plt.subplots(2,1,sharex=true,figsize=[3,5])
plt.subplots_adjust(wspace=0, hspace=0)
for i=1:dim
    if i in Calculation["method"].states
        axs[1,1].plot(r,Diabatic_Objects["potential"][:,i,i],label="V"*string(i))
        axs[1,1].plot(r,Objects["potential"][:,i,i],"--")
    end
end

axs[2,1].plot(r,Objects["nac"][:,1,2])
#

# axs[1,1].plot(r,-Diabatic_Objects["lx"][:,16,14],label="dia 13|"*string(16))
# axs[1,1].plot(r,-Objects["lx"][:,16,14],"--",label="adi 13|"*string(16))
# #
# axs[1,1].plot(r,-Diabatic_Objects["lx"][:,17,14],label="dia 13|"*string(17))
# axs[1,1].plot(r,-Objects["lx"][:,17,14],"--",label="adi 13|"*string(17))

# #
# axs[1,1].legend()
# #
# for i=1:dim
#     for j=i+1:dim
#         if (i in Calculation["method"].states)&(j in Calculation["method"].states)
#             axs[2,1].plot(r,Objects["nac"][:,i,j])
#             # if (lowercase(Calculation["method"].diabatisation) == "evolution")
#             #     axs[2,1].plot(r,[Objects["regularised_nac"][idx][i,j] for idx=1:lastindex(r)],"--",label="reg $i $j")
#             # end
#         end
#     end
# end
# axs[2,1].set_xlabel("Bond Length, Angstrom")
# axs[2,1].set_ylabel("NAC, 1/Angstrom")
# axs[1,1].set_ylabel("Potential, cm-1")
# axs[2,1].legend()


# print(Objects["potential"][1,:,:])


# save_diabatisation(Objects, Diabatic_Objects, lowercase(Calculation["method"].diabatisation), input_properties, "O2_3Pi", special_name = "NACs ")

# plt.figure()
# plt.plot(r,Objects["spin-orbit"][:,1,3],"--")
# plt.plot(r,Diabatic_Objects["spin-orbit"][:,1,3])

# plt.figure()
# plt.plot(r,Objects["lx"][:,1,3],"--")
# plt.plot(r,Diabatic_Objects["lx"][:,1,3])

# plt.figure()
# plt.plot(r,Objects["dipole"][:,1,3],"--")
# plt.plot(r,Diabatic_Objects["lx"][:,1,3])

# plt.plot(r,Objects["spin-orbit"][:,2,3],"--")
# plt.plot(r,Diabatic_Objects["spin-orbit"][:,2,3])

# plt.xlabel("Bond Length, Angstroms")
# plt.ylabel("Diabatic (potential) coupling, cm-1")