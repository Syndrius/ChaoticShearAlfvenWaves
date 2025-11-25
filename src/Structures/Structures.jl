"""
This module stores some of the key structures used throughout the package.
"""
module Structures


using FunctionWrappers
import FunctionWrappers: FunctionWrapper
using BSplineKit 


#good
include("Geometry.jl")

export GeometryT, MetT

#good
include("Island.jl") 

export IslandT, RadialIslandT, FluxIslandT, CoordIslandT
export init_island
export no_isl


#good
include("Fields.jl")

export FieldsT, BFieldT, FluxFieldsT, RadialFieldsT, IslandFieldsT

#good
include("FLR.jl")

export FLRT
export ideal_flr, init_flr


#good
include("Problem.jl")

export ProblemT

#good
include("QFM.jl")

export SurfaceITPT, TempSurfT, QFMSurfaceT

#good
include("Output.jl")

export EvalsT, find_ind

end
