"""

Module for the geometry. This includes
 - The main metric struct and functions for computing the different metric tensors.
 - The struct storing the island information.

"""

module Geometry


include("Metric.jl")

export MetT
export toroidal_metric!
export no_delta_metric!
export diagonal_toroidal_metric!
export cylindrical_metric!


include("Island.jl")

export IslandT


end