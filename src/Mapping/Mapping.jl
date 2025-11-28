"""

Module for mapping solution obtained in one coordinate system into another coordinate system for comparisons.
"""
module Mapping


using ..Structures
using ..Geometry 
using ..Grids
using ..QFM 
using ..PostProcessing
using ..Io
using ..Fields


using Elliptic
using JLD2
using Printf
using NLsolve 


#good enough
include("CoordTransform.jl")

#good enough.
include("Eigenfunctions.jl")

#good, and has been tested.
include("QFMToTor.jl")

export qfm_spectrum_to_tor

include("QFMToIsl.jl")

export qfm_spectrum_to_isl

#good and tested.
include("TorToIsl.jl")

#good and tested.
export tor_spectrum_to_isl

#good, and tested.
include("IslToTor.jl")

export isl_spectrum_to_tor

end
