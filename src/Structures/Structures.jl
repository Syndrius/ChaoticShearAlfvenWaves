"""
Module storing some of the key structures used throughout the package.
"""
module Structures


import FunctionWrappers: FunctionWrapper
import BSplineKit: SplineExtrapolations.SplineExtrapolation


include("Geometry.jl")

export GeometryT, MetT


include("Island.jl") 

export IslandT, RadialIslandT, FluxIslandT, CoordIslandT
export init_island, no_isl


include("Fields.jl")

export FieldsT, BFieldT, FluxFieldsT, RadialFieldsT, IslandFieldsT


include("FLR.jl")

export FLRT
export ideal_flr, init_flr


include("Problem.jl")

export ProblemT


include("QFM.jl")

export SurfaceITPT, TempSurfT, QFMSurfaceT


include("Output.jl")

export EvalsT, find_ind


end
