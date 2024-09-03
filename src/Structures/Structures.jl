"""

Contains the structures and constructors used, includes
 - The ProblemT type.
 - The GridsT type.
 - The EvalsT type.

"""
module Structures


using MID.Geometry
using MID.MagneticField



include("Grids.jl")

export GridsT
export FSSGridsT
export FFSGridsT
export FFFGridsT
export ContGridsT
export FEMGridDataT
export SMGridDataT
export mode_label
export init_fem_grid
export init_sm_grid
export init_grids
export instantiate_grids
export compute_ifft_grid


include("Problem.jl")

export ProblemT
export GeoParamsT
export FLRT
export init_problem


include("Output.jl")

export EvalsT
export find_ind


end