"""

Deals with the Hermite basis functions in 1, 2 and 3 dimensions for the finite element method. 
"""
module Basis


using ..Structures


include("ShapeFunctions.jl")



include("Hermite.jl")

export HB1dT, HB2dT, HB3dT
export hermite_basis, grid_id, basis_id


include("TrialFunction.jl")

export update_trial_function!


end
