"""

Module for mapping solution obtained in one coordinate system into another coordinate system for comparisons.
"""
module Mapping


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

#this module will require MIDParallel to be properly tested.

#good enough
include("CoordTransform.jl")


#good enough.
include("Eigenfunctions.jl")


#going to be a real prick to make sure these all work.
#probably need to check with the parallel version for all of these to work properly!

#export qfm_spectrum_to_tor

include("QFMToTor.jl")

export qfm_spectrum_to_tor

include("QFMToIsl.jl")

export qfm_spectrum_to_isl

include("TorToIsl.jl")

export tor_spectrum_to_isl

include("IslToTor.jl")

export isl_spectrum_to_tor

end
