"""
Module for solving the generalised eigenvalue problem PΦ = ω^2 QΦ.
"""
module Solve


using SparseArrays
using LinearAlgebra
using Arpack

using ..Structures


abstract type SolverT end


export solve, SolverT, init_solver

#good
include("Init.jl") 

#good
include("FullSpectrumSolve.jl")

#good
include("SliceSolve.jl")

#good
include("ShiftInvertSolve.jl")

end
