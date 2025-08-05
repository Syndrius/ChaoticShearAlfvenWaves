"""

Deals with the Hermite basis functions in 1, 2 and 3 dimensions and modifies the trial and test functions based on the local derivatives.
"""
module Basis

#Maybe consider moving some of MatrixSize and LocalBasis to Indexing. -> probably only if we ever generalise this code to allow different basis functions,

using ..Structures



include("MatrixSize.jl")

export local_matrix_size 
export matrix_size 
export init_basis_function


include("ShapeFunctions.jl")


include("Hermite.jl")

export hermite_basis



include("LocalBasis.jl")

export create_global_basis!
export local_to_global #arguably belongs in indexing...



end
