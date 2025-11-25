"""

Deals with the Hermite basis functions in 1, 2 and 3 dimensions for the finite element method. 
"""
module Basis


using ..Structures


#good
include("ShapeFunctions.jl")

#good
include("Hermite.jl")

export hermite_basis, HB1dT, HB2dT, HB3dT, grid_id, basis_id

#good
include("TrialFunction.jl")

export update_trial_function!


end
