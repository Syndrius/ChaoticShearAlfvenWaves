
"""

Quadratic Flux Minimisation. Module for creating approximate flux surfaces in chaotic fields, creating coordinates based on these surfaces and transforming variables into these coordinates.
"""

#TODO
#wrapping field lines is still a mystery. Few other smaller things to clean up, but otherwise this module is ok now.
#part of this is the rfft functions used by wrap_field_lines
#going to change this to just interacting with the surfaces rather than actually creating them, 
#creation will be moved to MIDCantori.

module QFM

using ..Structures
using ..Geometry
using ..Equilibrium


#using NLsolve
#using Roots
#using FFTW
using BSplineKit
using LinearAlgebra
#using Printf


#include("Action.jl")
#include("GradAction.jl")
#include("Rfft.jl")



include("Surfaces.jl")

export convert_surf
export compute_jac
export create_surf_itp


include("Transform.jl")

export CoordTransformT
export coord_transform!
export B_transform!
export met_transform!


end
