"""

Module for miscellaneous functions, predominatly around grids, indexing and matrix structure.

"""
module Misc


include("GIM.jl") #grids, indexing and Matrix structure

export compute_boundary_inds
export create_local_basis!
export grid_to_index
export index_to_grid
export reconstruct_phi
export clustered_grid
export cont_grid_to_index
export compute_length


include("FiniteElements.jl")

export RDataT
export hermite_basis
export local_to_global
export gauss_integrate
export gauss_integrate_for_big


include("Spectral.jl")

export ModeDataT
export spectral_grid

end