"""

Module for indexing our matrices. Deals with the mapping between 3d space and the matrix. Split into the different cases for different grids.
"""
module Grids


#import ..ChaoticShearAlfvenWaves: GridsT

using ..Basis

abstract type GridsT end #do we want this...?
abstract type GridT end


include("ContGrid.jl") #fine, probably needs to be incorporated properly
include("SMGrid.jl") #mostly fine will change if we modify the SM basis.
include("FEMGrid.jl") #fineish
include("ThreeDGrids.jl")

export GridsT, FSSGridsT, FFSGridsT, FFFGridsT, ContGridsT
export FEMGridT, SMGridT
export init_grids, init_grid
export inst_grids, mode_list, mode_label, periodic_grid, ifft_size, inst_grid


include("Boundaries.jl") #either fine, or needs a complete rewrite, depending on desire for flexible boundaries.

export compute_boundary_inds


include("Indexing.jl") #same as boundaries, currently fine, but if we want flexibility, this will all need to be rewritten.

export grid_to_index
export index_to_grid


#include("Basis.jl") #kind of weird, but will allow future changes more easily.

#export create_basis


include("MatrixSize.jl") #could be generalised if we take that approach

export matrix_size, init_trial_function, init_local_matrix

#this is a stupid name!
include("LocalBasis.jl")

export local_to_global!, local_to_global #not super sure this should be in here.

#probably a sign this module has gone to far tbh!
#perhaps a trial function module on its own? 
#this is always going to be awkward due to needing the basis, grids and arguably the weak form for this,
include("TrialFunction.jl")

export update_trial_function!


end
