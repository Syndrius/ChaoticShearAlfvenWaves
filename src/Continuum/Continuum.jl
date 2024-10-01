"""

Computes the continuum in the normal case with Br=0 and the specific island case with a linear rotational transform, as outlined by Qu and Hole 2022.
"""
module Continuum

using MID.Geometry
using MID.MagneticField
using MID.Structures
using MID.Indexing
using MID.WeakForm


using Elliptic 
using FFTW
using LinearAlgebra



struct psiIslandT
    m0 :: Int64 
    n0 :: Int64
    A :: Float64 
    q0 :: Float64 
    qp :: Float64 
    ψ0 :: Float64 
    w :: Float64 
end



include("ComputeContinuum.jl")

export continuum


include("IslandCoordinates.jl")

export compute_ψ̄


include("IslandWeakForm.jl")



include("IslandContinuum.jl")

export island_continuum

end