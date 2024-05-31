"""

Contains the structures of our inputs to our main functions, including
 - The ProblemT type.
 - The GridsT type.
 - File IO, for our inputs and for results.

"""
module Inputs


using MID.Geometry
using MID.MagneticField
using MID.Misc
using DelimitedFiles 


include("Grids.jl")

export GridsT
export init_grids
export construct_rgrid


include("Problem.jl")

export ProblemT
export GeoParamsT
export init_problem


include("FromFile.jl") #may want to move this to an Io module or something.

export inputs_to_file
export inputs_from_file
export eigvals_to_file 
export eigfuncs_to_file

end