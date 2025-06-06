"""

Module for the geometry. This includes
 - The main metric struct and functions for computing the different metric tensors.
 - The struct storing the island information.
 - Transformation between physical radius and toroidal flux.

"""
module Geometry

using Roots
using Elliptic
using LinearAlgebra



include("Island.jl")

export IslandT
export ContIslandT
export inst_island
export init_island
export sepratrix
export compute_sepratrix



include("Metric.jl")

export MetT
export toroidal_metric!
export no_delta_metric!
export diagonal_toroidal_metric!
export flux_toroidal_metric!
export cylindrical_metric!
export slab_metric!
export island_metric!



include("Flux.jl")

export f2r
export r2f


end
