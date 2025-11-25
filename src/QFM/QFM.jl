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

#good
include("Surfaces.jl")

export convert_surf
export compute_jac
export create_surf_itp


#good
include("Transform.jl")

export CoordTransformT
export coord_transform!
export B_transform!
export met_transform!


end
