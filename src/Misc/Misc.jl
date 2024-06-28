#TODO -> Need to think of a better name for this module. Currently handles basis indexing and numerical integration.
"""
This module will be rewritten :(

Module for miscellaneous functions, predominatly around grids, indexing and matrix structure.

"""
module Misc

using MID.Inputs #used for reconstructn=ing phi, not sure if it should belong here or not! Sort of makes sense

include("Basis.jl")

export hermite_basis
export create_local_basis!


include("IndexingFSS.jl") 
include("IndexingFFS.jl")

export compute_boundary_inds
export matrix_dim
export grid_to_index
export index_to_grid
export reconstruct_phi #should be the same eventually, not currently! -> shouldn't be here!
export cont_grid_to_index
export compute_length #not being used anymore tbh!
export local_to_global


include("Integration.jl")

export gauss_integrate

end