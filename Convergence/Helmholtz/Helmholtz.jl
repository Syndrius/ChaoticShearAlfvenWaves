
#testing module for solving the Helmholtz equation to ensure the finite elements is working as expected.
module Helmholtz

using FastGaussQuadrature
using LinearAlgebra
using ..Basis
using ..Structures
using ..Geometry
using ..Equilibrium
#using ..WeakForm
using Arpack


include("Structures.jl")
include("Basis.jl")
include("Indexing.jl")
include("OneD.jl")
include("TwoD.jl")
include("ThreeD.jl")
#include("WeakForm.jl")

export W_and_I!


end
