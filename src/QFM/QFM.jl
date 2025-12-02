"""

Module for converting to Quadratic Flux Minimising coordinates from QFM surfaces.
QFM surfaces are computed in the companion package CSAWCantori.
"""
module QFM


using ..Structures
using ..Grids
using ..Geometry 
using ..Fields


using BSplineKit
using LinearAlgebra


include("Surfaces.jl")

export convert_surf, compute_jac, create_surf_itp


include("Transform.jl")

export CoordTransformT
export coord_transform!, B_transform!, met_transform!


end
