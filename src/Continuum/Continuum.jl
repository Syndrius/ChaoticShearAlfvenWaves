"""

Computes the continuum in the normal case with Br=0 and the specific island case with a linear rotational transform, as outlined by Qu and Hole 2022.
"""
module Continuum

using MID.Geometry
using MID.MagneticField
using MID.Inputs
using MID.Indexing
using MID.WeakForm


using Elliptic 
using FFTW
using LinearAlgebra


include("ComputeContinuum.jl")

export continuum


include("IslandCoordinates.jl")

export compute_ψ̄


include("IslandWeakForm.jl")



include("IslandContinuum.jl")

export island_continuum

end