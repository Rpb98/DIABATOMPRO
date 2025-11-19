module DIABATOMPRO
#
## add PyPlot library, but be warned this is broken in some julia versions
using PyPlot
const plt = PyPlot
#
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
export Diabatise, fit_abinitio, Forward_Evolution, Backward_Evolution
#
function Diabatise(fname; save_flag = false, special_name = "")
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
    read_file(fname)
    #~~~~~~~~~~~~~~~~~~~~~~~~~ RUN HAMILTONIAN BUILDER ~~~~~~~~~~~~~~~~~~~~~~~~#
    include(joinpath(@__DIR__, "Build_Hamiltonian_Matrix.jl"))
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
    if Calculation["method"].abinitio_fit == true
        fit_abinitio()
        return nothing, nothing, nothing, nothing, nothing
    elseif Calculation["method"].abinitio_fit == false
        U, dU, UdU, K_Matrix, diabatic_basis, Diabatic_Objects, input_properties = run_diabatiser(lowercase(Calculation["method"].diabatisation))
    end
    #
    ## make a new Hamiltonian with different representations
    #
    ## save diabatisation?
    if save_flag == true
        save_diabatisation(Objects, Diabatic_Objects, lowercase(Calculation["method"].diabatisation), input_properties, fname, special_name = special_name)
    end
    #
    ## return the AtDT, adiabatic and diabatic matrix objects, diabatic basis, and Hamiltonians
    return U, Objects, Diabatic_Objects, diabatic_basis, Hamiltonian
end
end





