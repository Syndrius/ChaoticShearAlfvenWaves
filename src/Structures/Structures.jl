
module Structures


using FunctionWrappers
import FunctionWrappers: FunctionWrapper
using BSplineKit #not sure if this belongs!


#think this should be split into 3 files:
#WeakForm -> all structs needed for the weakform -> including Problem!
#Solver -> structs needed for solving!
#Grids
#Fields.

include("Geometry.jl")

export GeometryT, MetT, init_geometry

#unsure where this belongs tbh!
include("Island.jl") #TODO BIGTIME

export IslandT, RadIslandT, FluxIslandT, CoordIslandT

include("Fields.jl")

export FieldsT, BFieldT, init_fields, FluxFieldsT, RadFieldsT

include("FLR.jl")

export FLRT

include("Problem.jl")

export ProblemT, init_problem

include("Inputs.jl")

export WeakFormInputsT

include("QFM.jl")

export SurfaceITPT, TempSurfT, QFMSurfaceT

end
