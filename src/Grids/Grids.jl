"""

Module for creating and working with the grids.
Each dimension is defined individually then combined into a 3d grid.
"""
module Grids


using ..Basis

#abstract 1d grid.
abstract type GridT end


include("ContGrid.jl") 
include("SMGrid.jl") 
include("FEMGrid.jl") 
#good
include("ThreeDGrids.jl")

export GridsT, FSSGridsT, FFSGridsT, FFFGridsT, ContGridsT, FEMGridT, SMGridT
export init_grids, init_grid, inst_grid, inst_grids, mode_list, mode_label, periodic_grid, ifft_size


include("Boundaries.jl") 

export compute_boundary_inds


include("Indexing.jl") 

export grid_to_index, index_to_grid


include("MatrixSize.jl") 

export matrix_size, init_trial_function, init_local_matrix


include("LocalBasis.jl")

export local_to_global!


include("Interpolation.jl")

export interpolation


end
