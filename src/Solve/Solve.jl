"""

Module for solving the generalised eigenvalue problem PΦ = ω^2 QΦ.
Solving can be done
 - directly via Julia's LinearAlgebra, which is slow and not practical for large grids.
 - using shift and invert to target a specific frequency, obtaining the nev::Int64 nearest eigenvalues.
 - using a 'slicing' method where the shift invert method is used multiple times to build up a larger portion of the spectrum.
"""
module Solve


using SparseArrays
using LinearAlgebra
using Arpack


using ..Structures


abstract type SolverT end


export SolverT


include("Init.jl") 

export init_solver


include("FullSpectrumSolve.jl")
include("SliceSolve.jl")
include("ShiftInvertSolve.jl")

export ShiftInvertSolverT, FullSpectrumSolverT, SliceSolverT
export solve

end
