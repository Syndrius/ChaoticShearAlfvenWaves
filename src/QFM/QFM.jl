
"""

Quadratic Flux Minimisation. Module for creating approximate flux surfaces in chaotic fields, creating coordinates based on these surfaces and transforming variables into these coordinates.
"""

#this module is in a semi working state now
#TODO
#Clarity around Action.jl is terrible, need to understand the alg.
#Action grad is better but still pretty rough. Still allocating a lot.
#Rfft is not good.


module QFM

using ..Structures
using ..Geometry
using ..Equilibrium


using NLsolve
using FFTW
using BSplineKit
using LinearAlgebra
using Printf


include("Action.jl")
include("Rfft.jl")



include("Surfaces.jl")

export QFMSurfaceT
export SurfaceITPT
export TempSurfT #terrible name
export construct_surfaces
export farey_tree
export convert_surf
export create_surf_itp


include("Transform.jl")

export CoordTsfmT
export coord_transform!
export B_transform!
export met_transform!


end
