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
export RadIslandT
export FluxIslandT
export CoordIslandT
export ContIslandT
export inst_island
export init_island
export sepratrix
export compute_sepratrix
export convert_isl


#TODO
include("Metric.jl")
include("IslandMetric.jl")

export MetT
export metric!
export rad_toroidal_metric!
export flux_toroidal_metric!
export no_delta_metric!
export diagonal_toroidal_metric!
export flux_toroidal_metric!
export rad_cylindrical_metric!
export flux_cylindrical_metric!
export slab_metric!
export island_metric!



include("FluxConversion.jl")

export f2r
export r2f
export convert_island


end
