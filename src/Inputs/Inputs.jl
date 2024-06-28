"""

Contains the structures and constructors of our inputs to our main functions, including
 - The ProblemT type.
 - The GridsT type.

"""
module Inputs


using MID.Geometry
using MID.MagneticField


include("Grids.jl")

export GridsT
export FSSGridsT
export FFSGridsT
export FEMGridDataT
export SMGridDataT
export init_fem_grid
export init_sm_grid
export init_grids
export instantiate_grids


include("Problem.jl")

export ProblemT
export GeoParamsT
export init_problem


end