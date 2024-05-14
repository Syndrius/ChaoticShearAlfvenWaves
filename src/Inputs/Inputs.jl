

#defines the inputs needed to solve the problems.
module Inputs

using DelimitedFiles #makes me think this should be somewhere else.

using MID.Geometry
using MID.MagneticField
using MID.Misc


include("Grids.jl")

export GridsT
export init_grids


include("Problem.jl")

export ProblemT
export GeoParamsT
export init_problem


include("FromFile.jl") #may want to move this to an Io module or something.

export inputs_to_file
export inputs_from_file
export eigvals_to_file #may require reading these bad bois as well!
export eigfuncs_to_file

end