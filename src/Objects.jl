#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DEFINE OBJECT CLASSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
## create the PEC class
struct PEC
    ID           :: Int
    name         :: String
    obj_type     :: String
    mult         :: Number
    lambda       :: Number
    symmetry     :: String
    type         :: String
    sub_type     :: Tuple{String,String,String}
    units        :: Tuple{String,String}
    nparameters  :: Tuple{Int,Int,Int}
    Lval         :: Array{Any}
    Rval         :: Array{Number}
    factor       :: Float64
end
#
## create the spin-orbit class
struct SOC
    ID       :: Array{Int}
    name     :: String
    obj_type :: String
    spin     :: Tuple{Number,Number}
    sigma    :: Tuple{Number,Number}
    lambda   :: Tuple{Number,Number}
    Lz       :: Tuple{String,String}
    type     :: String
    units    :: Tuple{String,String}
    factor   :: Any
    Lval     :: Array{Any}
    Rval     :: Array{Number}
    sub_type :: Tuple{String,String,String}
end
#
## create the electronic angular momentum class
struct LX
    ID       :: Array{Int}
    name     :: String
    obj_type :: String
    spin     :: Tuple{Number, Number}
    lambda   :: Tuple{Number, Number}
    Lz       :: Tuple{String, String}
    type     :: String
    units    :: Tuple{String, String}
    factor   :: Any
    Lval     :: Array{Any}
    Rval     :: Array{Number}
    sub_type :: Tuple{String,String,String}
end 
#
## create the dipole moment class
struct DM
    ID       :: Array{Int}
    name     :: String
    obj_type :: String
    spin     :: Tuple{Number, Number}
    lambda   :: Tuple{Number, Number}
    Lz       :: Tuple{String, String}
    type     :: String
    units    :: Tuple{String, String}
    factor   :: Any
    Lval     :: Array{Any}
    Rval     :: Array{Number}
    sub_type :: Tuple{String,String,String}
end
#
## Create the NAC class
mutable struct NAC
    ID                :: Array{Int}
    name              :: String
    obj_type          :: String
    spin              :: Tuple{Number,Number}
    lambda            :: Tuple{Number,Number}
    type              :: String
    sub_type          :: Tuple{String,String,String}
    units             :: Tuple{String,String}
    factor            :: Any
    Lval              :: Array{Any}
    Rval              :: Array{Any}
    fit               :: Array{Int}
    bounds            :: Array{Any}
    fitted_parameters :: Array{Any}
    fit_range         :: Tuple{Number,Number}
end
#
## Create a switching function class
mutable struct SWITCH
    ID                :: Array{Int}
    obj_type          :: String
    type              :: String
    sub_type          :: Tuple{String,String,String}
    Lval              :: Array{Any}
    Rval              :: Array{Any}
    fit               :: Array{Int}
    bounds            :: Array{Any}
    fitted_parameters :: Array{Any}
    fit_range         :: Tuple{Number,Number}
end
#
## Create a grid class
struct Grid
    npoints            :: Int
    range              :: Tuple{Number,Number}
    interpolation_type :: String
end
#
## Create a method (of diabatisation) class
struct Method
    atoms                :: Array{Any}
    states               :: Array{Int}
    min_peak_distance    :: Number
    min_prominence       :: Float64
    thresh               :: Number
    diabatisation        :: String
    grid_resolution      :: Float64
    l_boundary_condition :: Any
    r_boundary_condition :: Any
    regularisation       :: Any
    plot                 :: Bool
    shift                :: Bool
    abinitio_fit         :: Bool
end
#
## create a diabatisation save class
struct Save
	as   :: String
	jrot :: String
	vmax :: String
	nroots :: String
end
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ OBJECT CLASS METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
## define operators
Base.isless(a::PEC, b::PEC) = a.ID < b.ID
Base.isgreater(a::PEC, b::PEC) = a.ID > b.ID