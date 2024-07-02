"""

Deals with the Hermite basis functions in 1, 2 and 3 dimensions and modifies the trial and test functions based on the local derivatives.
"""
module Basis


using MID.Inputs


include("FSSBasis.jl")
include("FFSBasis.jl")
include("FFFBasis.jl")

export hermite_basis 
export create_local_basis!


end