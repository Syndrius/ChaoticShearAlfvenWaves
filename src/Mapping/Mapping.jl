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



include("CoordTransform.jl")


include("Eigenfunctions.jl")


include("QFMToTor.jl")

export qfm_spectrum_to_tor


include("QFMToIsl.jl")

export qfm_spectrum_to_isl


include("TorToIsl.jl")

export tor_spectrum_to_isl


include("IslToTor.jl")

export isl_spectrum_to_tor

end
