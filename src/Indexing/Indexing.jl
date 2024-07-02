"""

Module for indexing our matrices. Deals with the mapping between 3d space and the matrix. Split into the different cases for different grids.
"""
module Indexing


using MID.Inputs


"""
Arrays that distinguish the different Hermite basis functions, used for converting between the grid and the appropriate points indices in the matrices.
"""
const grid_id = [0, 0, 1, 1]
const basis_id = [0, 1, 0, 1]


include("FSSIndexing.jl")
include("FFSIndexing.jl")
include("FFFIndexing.jl")


export compute_boundary_inds
export matrix_dim
export grid_to_index
export index_to_grid
export reconstruct_phi 
export cont_grid_to_index
export local_to_global


end