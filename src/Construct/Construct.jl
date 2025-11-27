"""

Module for Constructing the two main matrices in the equation PΦ = ω²QΦ. 
Construction is split up based on which grids are used. 
Additionally, construction with qfm surfaces is slightly modified.
"""
module Construct

using FFTW
using SparseArrays
using LinearAlgebra


using ..Geometry
using ..Basis 
using ..Grids
using ..Integration
using ..Fields
using ..WeakForm
using ..Structures
using ..QFM


include("FSS.jl")
include("FFS.jl")
include("FFF.jl")

export construct


include("Continuum.jl")

export continuum



end
