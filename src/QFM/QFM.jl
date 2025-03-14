
"""

Quadratic Flux Minimisation. Module for creating approximate flux surfaces in chaotic fields, creating coordinates based on these surfaces and transforming variables into these coordinates.
"""

#this module is in a semi working state now
#TODO
#determine example that matches better with our setup
#add bounding surfaces and I guess the straightening
#better integration into MID properly
#speed up finding the flux surface/fix the functions a bit
#perhaps add the derivative of the action or whatever it is for I guess better solving?
#probably see if we can do some example of this with minimally modified flux surfaces to see if our results are ok then.

module QFM

using ..Structures
using ..Geometry
using ..MagneticField

#using MID.Indexing
#using MID.Basis


using NLsolve
using FFTW
#using LaTeXStrings #note only here for plotting which should be removed eventually!
#using Plots #same with this!!!
#using Dierckx #don't think this will work
#using Interpolations
using BSplineKit
using LinearAlgebra
#unsure where this will go in the heirachy yet!

#this needs some serious work still
#action1 is probably obsolete, action2 is better, but grad_Jm is a (working!) disaster, needs to be cleaned and optimised
include("Action.jl")
include("Action2.jl")

#don't think this needs to export anything! could be a mistake though!


include("Surfaces.jl")

export QFMSurfaceT
export SurfaceITPT
export construct_surfaces
export farey_tree
export convert_surf
export create_surf_itp

#may need to split this up.
include("Transform.jl")

export CoordTsfmT
export coord_transform!
export B_transform!
export met_transform!




end