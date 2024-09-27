"""

Deals with the Hermite basis functions in 1, 2 and 3 dimensions and modifies the trial and test functions based on the local derivatives.
"""
module Basis


using MID.Structures



include("MatrixSize.jl")

export local_matrix_size 
export matrix_size 
export init_bases_function



include("Hermite.jl")

export hermite_basis



include("LocalBasis.jl")

export create_local_basis!
export local_to_global #arguably belongs in indexing...




end