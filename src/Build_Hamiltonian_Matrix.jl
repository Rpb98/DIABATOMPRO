#~~~~~~~~~~~~~~~~~~~~~~ BUILDING THE HAMILTONIAN MATRIX ~~~~~~~~~~~~~~~~~~~~~~~#
#
## order the Potential Dictionary on the ID's
# Potential = Dict(sort(collect(ioRead.Potential)))
function build_hamiltonian_objects(Calculation, Potential, SpinOrbit, EAMC, Dipole, NonAdiabaticCoupling)
    #
    ## define dimension of system (how many PECs)
    dim = length(Potential)
    #
    ## intitialise a matrix to populate the numerical potentials into by column
    PotMat = zeros(Calculation["grid"].npoints,dim,dim) #Array{Float64}(undef, Calculation["grid"].npoints,dim,dim) 
    SOMat  = zeros(Calculation["grid"].npoints,dim,dim) #Array{Float64}(undef, Calculation["grid"].npoints,dim,dim) 
    DipMat = zeros(Calculation["grid"].npoints,dim,dim) #Array{Float64}(undef, Calculation["grid"].npoints,dim,dim)
    EAMMat = zeros(Calculation["grid"].npoints,dim,dim) #Array{Float64}(undef, Calculation["grid"].npoints,dim,dim)
    NACMat = zeros(Calculation["grid"].npoints,dim,dim) #Array{Float64}(undef, Calculation["grid"].npoints,dim,dim)
    #
    ##
    X = Potential[1].Rval
    Xre  = minimum(X)
    # print(Xre)
    #
    ## loop over potential instances and compute the numerical potentials
    for key in keys(Potential)
        # if any(Calculation["method"].states.==key)
        local r, V
        r, V = ComputeProperty(Potential[key])
        #
        PotMat[:,key,key] = V #.- Xre*219474.6313708000
        # end
    end
    #
    ## SOCs
    for key in keys(SpinOrbit)
        local r, SO
        r, SO = ComputeProperty(SpinOrbit[key])
        SOMat[:,key[1],key[2]] =  SO
        #
        ## now set the dagger element in the matrix if not specified in input and
        ## not a diagonal term
        if (key[1] != key[2])&([key[2],key[1]] ∉ keys(SpinOrbit))
            SOMat[:,key[2],key[1]] = -SO
        end
    end
    #
    ## EAMCs
    for key in keys(EAMC)
        local r, EAM
        r, EAM = ComputeProperty(EAMC[key])
        EAMMat[:,key[1],key[2]] = EAM
        #
        ## now set the dagger element in the matrix if not specified in input and
        ## not a diagonal term
        if (key[1] != key[2])&([key[2],key[1]] ∉ keys(EAMC))
            EAMMat[:,key[2],key[1]] = -EAM
        end
    end
    #
    ## dipoles
    for key in keys(Dipole)
        local r, DM
        r, DM = ComputeProperty(Dipole[key])
        DipMat[:,key[1],key[2]] = DM
        #
        ## now set the dagger element in the matrix if not specified in input and
        ## not a diagonal term
        if (key[1] != key[2])&([key[2],key[1]] ∉ keys(Dipole))
            DipMat[:,key[2],key[1]] = -DM
        end
    end
    #
    ## NACs
    for key in keys(NonAdiabaticCoupling)
        local r, NAC
        r, NAC = ComputeProperty(NonAdiabaticCoupling[key])
        NACMat[:,key[1],key[2]] = NAC
        #
        ## now set the dagger element in the matrix if not specified in input and
        ## not a diagonal term
        if (key[1] != key[2])&([key[2],key[1]] ∉ keys(NonAdiabaticCoupling))
            NACMat[:,key[2],key[1]] = -NAC
        end
    end
    #
    ## make a dictionary of these matrices
    Objects = Dict()
    Objects["potential"]=PotMat
    Objects["spin-orbit"]=SOMat
    Objects["lx"]=EAMMat
    Objects["dipole"]=DipMat
    Objects["nac"]=NACMat
    #
    return Objects
end