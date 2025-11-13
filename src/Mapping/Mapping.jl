"""

Module for mapping solution obtained in one coordinate system into another coordinate system for comparisons.
"""
module Mapping

#We may one day want to move the hermite interpolation into Post-processing
#in theory we may want to plot with a finer grain or whatever using the hermite interpolation.

#think it is only using all of this shite because of struct definitions.

using ..Structures
using ..Geometry #probably don't need this? Perhaps for coords.
using ..Grids
using ..QFM 
using ..PostProcessing
using ..Io
using ..Basis #this is either fine, or basis should do the interpolation.


using Elliptic
using JLD2
using Printf
using NLsolve 


include("CoordTransform.jl")


include("HermiteInterpolation.jl") 


include("DerivativeInterpolation.jl") 


include("Eigenfunctions.jl")


include("Spectrum.jl") 


#probably just remove this.
include("Harmonics.jl") 



end
