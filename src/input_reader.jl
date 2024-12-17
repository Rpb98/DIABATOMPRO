#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INPUT FILE READER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
function display_progress(text,current, total)
    percent_complete = round(current / total * 100)
    bar_length = 29
    num_bar_symbols = Int(round(current / total * bar_length))  # Convert to integer
    bar = "" * repeat("█" , num_bar_symbols) * repeat(" ", bar_length - num_bar_symbols) * ""
    
    # Use carriage return to overwrite the previous line
    print("\r"*text*": $percent_complete% $bar ($current/$total)")
end
#
## comment remover
function remove_comment(str::AbstractString)::String
    # Define regular expression pattern to match substrings enclosed in 
    # parentheses
    pattern = r"\([^)]*\)"
    # Replace all matched substrings with an empty string
    return replace(str, pattern => "")
end
#
function message(str::AbstractString)
    #
    if occursin("morphing",str)
        # print(str,"\n")
        # print(" ***Morphing is not implemented as of version 1.0. \n")
        return true
    end
end
#
## Initialise the Hamiltonian and sub-Hamiltonians as a dictionary of object 
## classes
Hamiltonian          = Dict{Tuple{String, Any}, Any}()
Potential            = Dict{Int, Any}()
SpinOrbit            = Dict{Any, Any}()
EAMC                 = Dict{Any, Any}()
Dipole               = Dict{Any, Any}()
NonAdiabaticCoupling = Dict{Any, Any}()
SwitchingFunction    = Dict{Any, Any}()
#
## Intitialise the Calculation dictionary which defines hyper parameters
Calculation          = Dict{String, Any}()
#
## ab initio classes
abinitio = Dict{Any, Any}()
#
## make the Hamiltonians global quantities
global Hamiltonian, Potential, SpinOrbit, EAMC, Dipole, NonAdiabaticCoupling, Calculation, abinitio, SwitchingFunction
#
function process_line(line::String)::String
    #
    ## remove comments from line
    line = remove_comment(line)
    #
    ## Strip leading or trailing whitespace from line
    line = strip(line)
    #
    return line
end
#
function is_start_of_block(line::String, objects, read_flag::Bool)
    #
    ## define the beginning of the object block conditions
    #  - is the line only one object long?
    #  - does the line contain any of the objects (method, potential, NAC, e.t.c.)?
    #  - is the field not an abinitio field?
    #  - is the code currently reading the object field?
    #
    if  (size(split(line))[1]!=0)
        return (any(occursin.(objects,lowercase(split(line)[1]))))&&
            (read_flag==false)#&(occursin("abinitio",lowercase(line))==false)
    else
        return false
    end
end
#
function initialise_input_parameters(line::String)
    #
    ## Define a dictionary to hold the input values and initialise some default parameters
    input_values = Dict()
    #
    input_values["factor"] = "1.0"
    #
    input_values["units"]  = ("angstrom", "au")
    #
    input_values["symmetry"]  = "+"
    #
    input_values["<x|lz|y>"] = ("N/A","N/A")
    #
    input_values["sub-types"] = ("N/A","N/A","N/A")
    #
    input_values["nparameters"] = ("-1","-1","-1")
    #
    input_values["rref"] = "-1"
    #
    input_values["fit"] = []
    #
    input_values["bounds"] = []
    #
    input_values["shift"] = "true"
    #
    input_values["diabatisation"] = "evolution"
    #
    input_values["regularisation"] = false
    #
    input_values["l_boundary_condition"] = false
    #
    input_values["r_boundary_condition"] = false
    #
    input_values["grid_resolution"] = 0.0001
    #
    input_values["save"] = "false" 
    #
    input_values["plot"] = "false"
    #
    input_values["abinitio_fit"] = "false"
    #
    input_values["fit_range"] = ("-1e100","1e100")
    #
    ## save block
    input_values["as"] = "file"
    #
    input_values["jrot"] = "0"
    #
    input_values["vmax"] = "1"
    #
    input_values["nroots"] = "1"
    #
    ## legacy parameters
    input_values["min_peak_distance"] = "0.2"
    #
    input_values["min_prominence"] = "8e10"
    #
    input_values["thresh"] = "0.1"
    #
    ## process the current line for ID and object type
    #
    ## extract the subfields in the line: object type, IDs
    field = split(line)
    #
    object_key = field[1]
    ID = field[2:end]
    #
    ## Strip whitespace from the key and ID
    object_key = strip(object_key)
    ID = strip.(ID)
    #
    ## add the object to the input dictionary
    input_values["id"] = ID
    #
    return input_values, object_key, ID
end
#
function add_key_value_pair(input_values::Dict{Any,Any}, key::String, value::Any)
    #
    ## 
    # if (size(value)[1]>2)&(key!="states")&(key!="nparameters")&(key!="sub-types")&(key!="r_boundary_condition")&(key!="l_boundary_condition")
    #     print("cats ",key,"\n")
    #     #print("*** ERROR: too many parameters in the the line above ^^^. Maximum 2. \n")
    #     break
    #
    ## states involved in the object
    if key=="states"
        input_values[lowercase(key)] = value
    #
    ## sub-types, i.e. for functions, mixing angles, perturbation functions e.t.c.
    elseif key=="sub-types"
        st = ["N/A" for _=1:3]
        for (idx, val) in enumerate(value)
            st[idx] = val
        end
        input_values[lowercase(key)] = Tuple(st)
    #
    elseif key=="nparameters"
        input_values[lowercase(key)] = Tuple(value) #(value[1], value[2], value[3])
    #
    elseif occursin("boundary_condition",key)
        N = Int(sqrt(size(value)[1]))
        M = zeros(Float64,N,N)
        k=0
        #
        for i=1:N
            for j=1:N
                k+=1
                M[i,j] = parse(Float64,value[k])
            end
        end
        input_values[lowercase(key)] = M
    #
    elseif lowercase(key) == "jrot"
        input_values[lowercase(key)] = join(value)
    elseif lowercase(key) == "vmax"
        input_values[lowercase(key)] = join(value)
    #
    elseif size(value)[1] != 1
        input_values[lowercase(key)] = Tuple(value)
    #
    elseif size(value)[1] == 1
        input_values[lowercase(key)] = value[1]
    end
    #
    return input_values
end
#
function create_object_instance(input_values::Dict{Any,Any}, object_key)
    #
    ## grid instance
    if object_key == "grid"                                         
        #
        grid = Grid(parse(Int,input_values["npoints"]),
                          parse.(Float64,input_values["range"]),
                          input_values["interpolation_type"]
                          )
        #
        Calculation["grid"] = grid
    #
    ## method instance
    elseif object_key == "method"                                       
        #
        ## if no left boundary condition is given then use the identity matrix
        if input_values["l_boundary_condition"] == false
            Identity = zeros(Float64,length(input_values["states"]),length(input_values["states"]))
            for i=1:length(input_values["states"])
                Identity[i,i] = 1.0
            end
            #
            input_values["l_boundary_condition"] = Identity
        end
        #
        ## check if atoms have been given a numeric mass
        masses = [0.0,0.0]
        #
        for i=1:2
            mass_flag = tryparse(Float64, input_values["atoms"][i]) 
            if mass_flag !== nothing
                masses[i] = mass_flag
            else
                masses[i] = atomic_mass[input_values["atoms"][i]]
            end
        end
        #
        method = Method(masses,
                        parse.(Int,input_values["states"]),
                        parse.(Float64,input_values["min_peak_distance"]),
                        parse.(Float64,input_values["min_prominence"]),
                        parse.(Float64,input_values["thresh"]),
                        input_values["diabatisation"],
                        parse.(Float64,input_values["grid_resolution"]),
                        input_values["l_boundary_condition"],
                        input_values["r_boundary_condition"],
                        input_values["regularisation"],
                        parse.(Bool,input_values["plot"]),
                        parse.(Bool,input_values["shift"]),
                        parse.(Bool,input_values["abinitio_fit"])
                    )
        #
        Calculation["method"] = method
    #
    ## save diabatisation instance
    elseif object_key == "save"                                    
        #
        save = Save(input_values["as"],
                    input_values["jrot"],
                    input_values["vmax"],
                    input_values["nroots"]
                    )
        #
        Calculation["save"] = save
    #
    ## potential objects
    elseif object_key == "poten"                                  
        obj_type = "poten"
        #
        ## if object has grid type, parse left column to floats (bonds)
        if input_values["type"]=="grid"
            input_values["Lval"] = parse.(Float64,input_values["Lval"])
        end
        #
        ## 
        object = PEC(parse.(Int,input_values["id"])[1],
                                input_values["name"],
                                "PEC",
                 parse.(Float64,input_values["mult"]),
                 parse.(Float64,input_values["lambda"]),
                                input_values["symmetry"],
                                input_values["type"],
                                input_values["sub-types"],
                                input_values["units"],
                     parse.(Int,input_values["nparameters"]),
                                input_values["Lval"],
                                input_values["Rval"],
                  parse(Float64,input_values["factor"]))
        #
        Potential[object.ID]=object
        #
        return object, obj_type
    #
    ## spin-orbit objects
    elseif occursin("spin-orbit",lowercase(object_key))               
        obj_type = "spin-orbit"
        #
        ## if object has grid type, parse left column (r) to floats 
        if input_values["type"]=="grid"
            input_values["Lval"] = parse.(Float64,input_values["Lval"])
        end
        #
        ##
        object = SOC(parse.(Int,input_values["id"]),
                     input_values["name"],
                     "SOC",
                     parse.(Float64,input_values["spin"]),
                     parse.(Float64,input_values["sigma"]),
                     parse.(Float64,input_values["lambda"]),
                     input_values["<x|lz|y>"],
                     input_values["type"],
                     input_values["units"],
                     parse.(Float64,input_values["factor"]),
                     input_values["Lval"],
                     input_values["Rval"],
                     input_values["sub-types"])
        #
        SpinOrbit[object.ID]=object
        #
        return object, obj_type
    #
    ##
    elseif occursin("dipole",lowercase(object_key))                     # dipole instance
        obj_type = "dipole"
        #
        ## if object has grid type, parse left column to floats (bonds
        if input_values["type"]=="grid"
            input_values["Lval"] = parse.(Float64,input_values["Lval"])
        end
        #
        object = DM(parse.(Int,input_values["id"]),
                    input_values["name"],
                    "DM",
                    parse.(Float64,input_values["spin"]),
                    parse.(Float64,input_values["lambda"]),
                    input_values["<x|lz|y>"],
                    input_values["type"],
                    input_values["units"],
      parse(Float64,input_values["factor"]),
                    input_values["Lval"],
                    input_values["Rval"],
                    input_values["sub-types"]) 
        #
        Dipole[object.ID]=object
        #
        return object, obj_type
    #
    ##
    elseif object_key == "LX"                                           # dipole instance
        obj_type = "LX"
        #
        ## if object has grid type, parse left column to floats (bonds
        if input_values["type"]=="grid"
            input_values["Lval"] = parse.(Float64,input_values["Lval"])
        end
        #
        object = LX(parse.(Int,input_values["id"]),
                    input_values["name"],
                    "LX",
                    parse.(Float64,input_values["spin"]),
                    parse.(Float64,input_values["lambda"]),
                    input_values["<x|lz|y>"],
                    input_values["type"],
                    input_values["units"],
      parse(Float64,input_values["factor"]),
                    input_values["Lval"],
                    input_values["Rval"],
                    input_values["sub-types"]) 
        #
        EAMC[object.ID]=object
        #
        return object, obj_type
    #
    ## 
    elseif occursin("nac",lowercase(object_key))                                          # NAC instance
        obj_type = "NAC"
        #
        range_init = parse.(Float64,input_values["fit_range"])
        #
        if (range_init[1] == -1e100)&(range_init[2]==1e100)
            range = Calculation["grid"].range 
        else
        # if occursin("abinitio",lowercase(object_key))
            range = parse.(Float64,input_values["fit_range"])
        end
        #
        ## if object has grid type, parse left column to floats (bonds
        if input_values["type"]=="grid"
            input_values["Lval"] = parse.(Float64,input_values["Lval"])
        end
        #
        object = NAC(parse.(Int,input_values["id"]),
                     input_values["name"],
                     "NAC",
                     parse.(Float64,input_values["spin"]),
                     parse.(Float64,input_values["lambda"]),
                     input_values["type"],
                     input_values["sub-types"],
                     input_values["units"],
                     parse(Float64,input_values["factor"]),
                     input_values["Lval"],
                     input_values["Rval"],
                     input_values["fit"],
                     input_values["bounds"],
                     input_values["Rval"],
                     range) 
        #
        if occursin("abinitio",lowercase(object_key))
            abinitio[(obj_type,object.ID)]=object 
        else
            NonAdiabaticCoupling[object.ID]=object 
        end
        #
        return object, obj_type
    elseif lowercase(object_key) == "switch"
        obj_type = "switch"
        #
        object = SWITCH(parse.(Int,input_values["id"]),
                        "switch",
                        input_values["Lval"],
                        input_values["Rval"],
                        input_values["fit"],
                        input_values["bounds"],
                        input_values["Rval"]
                        ) 
        #
        SwitchingFunction[object.ID]=object
        # 
        return object, obj_type
    end
end
#
function values_column_pop(line::String, fitcol::Vector{Any}, boundscol::Vector{Any}, Lcol::Array{Any}, Rcol::Array{Any})
    #
    ## Split the line into a key and value on whitespace
    ln = split(line)
    Lval, Rval = ln[1], ln[2:end]
    Lval = strip(Lval)
    Rval = strip.(Rval)
    #
    ## depending on line length populate fit, bounds, L and R columns
    if (length(ln) >= 3)
        if (ln[3] == "fit")
            if (length(ln) == 3)
                push!(fitcol,1)
                push!(boundscol,[-1e20,1e20])
            elseif (length(ln) > 3)
                push!(fitcol,1)
                push!(boundscol,[parse(Float64,ln[4]),parse(Float64,ln[5])])
            end
            #
            Rval = parse(Float64,Rval[1])
        else
            push!(fitcol,0)
            push!(boundscol,[NaN,NaN])
            #
            Rval = parse.(Float64,Rval)
            if length(Rval) == 1
                Rval = Rval[1]
            end
        end
    else
        push!(fitcol,0)
        push!(boundscol,[NaN,NaN])
        Rval = parse.(Float64,Rval)
        if length(Rval) == 1
            Rval = Rval[1]
        end
    end
    #
    push!(Lcol,Lval)
    push!(Rcol,Rval)
end
#
## Open the input file for reading
function read_file(fname)
    open(fname) do f
        # numberOfLines = countlines(fname)
        #
        ## initialise parameters
        input_values = Dict()
        #
        ## flags:
        read_flag    = false     
        extract_flag = false 
        #
        ## object key words: (i.e. poten, spin-orbit,...)
        object_key = 0
        #
        ## Left and Right columns of object block values:
        Lcol = []
        Rcol = []
        fitcol = []
        boundscol = []
        #
        ## objects of input file
        objects = ["method",
                     "grid",
                     "save",
                    "poten",
               "spin-orbit",
                   "dipole", 
                       "lx", 
                      "nac",
                   "switch"]
        #
        ## Loop over each line in the input file
        for (LineIndex, line) in enumerate(eachline(f))
            #
            ## terminal print out messages with continue check
            if message(line) == true; continue; end
            #
            ## process line (whitespaces, comments, e.t.c.)
            line = process_line(line)
            #
            ## Initialise the object block extraction if at beginning of object field                   
            if is_start_of_block(line, objects, read_flag)                          # begin object extraction
                #
                ## initialise the input values dictionary with default values
                input_values, object_key, ID = initialise_input_parameters(line)        
                #
                ## switch the object read block flag on
                read_flag = true                                                    # read_flag = true
                continue
            end
            #
            ## begin reading the fields in the object block
            if read_flag == true                                                    # reading fields in object block
                #
                ## check if the end of the block has been reached
                if occursin("end",line)                                             # check 'end' of the block field
                    #
                    ## set read and exraction flags now we at the end of object
                    read_flag = false
                    extract_flag = false
                    #
                    ## populate the calculation dictionary object
                    create_object_instance(input_values, object_key)
                    continue
                end
                #
                ## check to see if the value sub-block has been reached
                if occursin("values",line)                                          # check to see if 'values' key in line
                    #
                    ## define the left and right column arrays defining the values
                    ## of the object
                    Lcol = []
                    Rcol = []
                    fitcol = []
                    boundscol = []
                    #
                    ## initialise value extraction flag and turn off read_flag
                    extract_flag = true                                             # extract_flag = true
                    read_flag = false                                               # read_flag = false
                    continue
                end
                #
                # Split the line into a key and value on whitespace
                field = split(line)
                key   = lowercase(field[1])
                value = field[2:end]
                #
                ## Add the key-value pair to the dictionary
                input_values = add_key_value_pair(input_values, key, value)
            end
            #
            if extract_flag ==  true                                                # begin extraction of object
                #
                ## first check to see if the end of the block has been reached
                if occursin("end",lowercase(line))                                  # end of object, create instance
                    #
                    ## set the object extraction flag to false
                    extract_flag = false
                    #
                    ## Add the bonds and unit vals to the dictionary
                    input_values["Lval"]    = Lcol
                    input_values["Rval"]    = Rcol
                    input_values["fit"]     = fitcol
                    input_values["bounds"]  = boundscol
                    #
                    ## Now create an instance of the object
                    object, obj_type = create_object_instance(input_values, object_key)
                    #
                    ## if the object block corresponded to a molecular property, 
                    ## insert its class instance into the Hamiltonian dictionary
                    if (object_key != "grid")&(object_key != "method")&(occursin("abinitio",lowercase(object_key)) == false)&(lowercase(object_key) != "switch") # populate Hamiltonian
                        Hamiltonian[(obj_type,object.ID)]=object
                    end
                    continue
                end
                #
                ## populate the Left, Right, fitting, and bounds column for grid/parameter values
                values_column_pop(line,fitcol, boundscol, Lcol, Rcol)
            end
        end
    end
end