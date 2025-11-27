"""

Module for creating and working with the grids.
This includes indexing between grids and matrices and creating the trial and test functions.
"""
module Grids


using ..Basis

#abstract 1d grid.
abstract type GridT end

#good
include("ContGrid.jl") 
#good
include("SMGrid.jl") 
#good
include("FEMGrid.jl") 
#good
include("ThreeDGrids.jl")

export GridsT, FSSGridsT, FFSGridsT, FFFGridsT, ContGridsT
export FEMGridT, SMGridT
export init_grids, init_grid, inst_grid, inst_grids
export mode_list, mode_label, periodic_grid, ifft_size

#good
include("Boundaries.jl") 

export compute_boundary_inds


#good
include("Indexing.jl") 

export grid_to_index
export index_to_grid

#good
include("MatrixSize.jl") 

export matrix_size, init_trial_function, init_local_matrix

#good
include("LocalBasis.jl")

export local_to_global!

include("Interpolation.jl")

export interpolation


end
