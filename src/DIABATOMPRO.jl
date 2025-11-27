module DIABATOMPRO
#
## add PyPlot library, but be warned this is broken in some julia versions
# using PyPlot
# const plt = PyPlot
#
struct DiabatisationResult
    U
    adiabatic
    diabatic
    input_properties
    residual_kin
    H
end

# Tell Julia how to print this object nicely
function Base.show(io::IO, x::DiabatisationResult)
    println(io, "✅ Diabatisation Complete.")
    println(io, "   • Results stored in output object. Access via:")
    println(io, "       .U = AtDT")
    println(io, "       .adiabatic = Adiabatic Property Dictionary")
    println(io, "       .diabatic = Diabatic Property Dictionary")
    println(io, "       .input_properties = list of Electronic strucure properties input by user")
    println(io, "       .residual_kin = vector of residual kinetic energy matrix norms")
    println(io, "       .H = Total Hamiltonian Dictionary")
end
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUN HAMILTONIAN BUILDER ~~~~~~~~~~~~~~~~~~~~~~~~~~#
# include(joinpath(@__DIR__, "Build_Hamiltonian_Matrix.jl"))
include("Build_Hamiltonian_Matrix.jl")
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
    print("𝚆𝚛𝚒𝚝𝚝𝚎𝚗 𝚋𝚢 𝚁𝚢𝚊𝚗 𝙱𝚛𝚊𝚍𝚢 (UCL) | 𝚟1.1 | 𝙻𝚊𝚜𝚝 𝚞𝚙𝚍𝚊𝚝𝚎𝚍 20.𝟷1.𝟸𝟶𝟸5\n")
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
    # 2. CALL THE FUNCTION HERE
    # We pass the global variables (Calculation, Potential, etc.) into the function
    # and capture the result in 'Objects'
    global Objects = build_hamiltonian_objects(Calculation, 
                                               Potential, 
                                               SpinOrbit, 
                                               EAMC, 
                                               Dipole, 
                                               NonAdiabaticCoupling)
    #
    # Warning for specific versions
    if VERSION < v"1.6.7"
        @warn "You are using Julia $VERSION. DIABATOMPRO is tested on 1.6.7, 1.12. Some features may fail."
    end
    #
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
        U, dU, UdU, K_Matrix, Diabatic_Objects, input_properties, residual_kinetic_energy = run_diabatiser(lowercase(Calculation["method"].diabatisation))
    end
    #
    ## make a new Hamiltonian with different representations
    #
    ## save diabatisation?
    if save_flag == true
        println("___SAVING____")
        save_diabatisation(Objects, Diabatic_Objects, lowercase(Calculation["method"].diabatisation), input_properties, fname, special_name = special_name)
    end
    #
    ## return the AtDT, adiabatic and diabatic matrix objects, diabatic basis, and Hamiltonians
    return DiabatisationResult(U, Objects, Diabatic_Objects, input_properties, residual_kinetic_energy, Hamiltonian)
end
end





