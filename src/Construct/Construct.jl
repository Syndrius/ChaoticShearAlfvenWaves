
module Construct

using FFTW
using SparseArrays
using FastGaussQuadrature
using LinearAlgebra


using ..Geometry
using ..Basis
using ..Indexing
using ..Integration
using ..Equilibrium
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
