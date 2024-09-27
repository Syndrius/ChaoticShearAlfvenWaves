"""

Module for indexing our matrices. Deals with the mapping between 3d space and the matrix. Split into the different cases for different grids.
"""
module Indexing


using MID.Structures


"""
Arrays that distinguish the different Hermite basis functions, used for converting between the grid and the appropriate indices in the matrices.
"""
const grid_id = [0, 0, 1, 1]
const basis_id = [0, 1, 0, 1]



include("Boundaries.jl")

export compute_boundary_inds



include("GridToIndex.jl")

export grid_to_index
export index_to_grid


#needs serious work, should conform to everything else.
#include("ContIndexing.jl")
#export cont_grid_to_index


end