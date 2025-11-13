
module Solve


#this is a much cleaner module, but we require some kind of normalisation in order to make this work.

using SparseArrays
using LinearAlgebra

using ..Structures

abstract type SolverT end


export solve, SolverT, init_solver

include("Init.jl") #probably want to move this into this file.

include("FullSpectrumSolve.jl")

include("SliceSolve.jl")


include("ShiftInvertSolve.jl")

end
