"""

Module for mapping solution obtained in one coordinate system into another coordinate system for comparisons.
"""
module Mapping

#We may one day want to move the hermite interpolation into Post-processing
#in theory we may want to plot with a finer grain or whatever using the hermite interpolation.

using ..Structures
using ..Geometry
using ..Indexing 
using ..QFM 
using ..PostProcessing
using ..Io


using Elliptic
using JLD2
using Printf
using NLsolve 


include("CoordTransform.jl")


include("HermiteInterpolation.jl") 


include("Eigenfunctions.jl")


include("Spectrum.jl") 



end
