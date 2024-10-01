"""

Contains the structures and constructors used, includes
 - The ProblemT type.
 - The GridsT type.
 - The EvalsT type.

"""
module Structures


using MID.Geometry
using MID.MagneticField


"""
Abstract type for each grid.
"""
abstract type GridDataT end


"""
Abstract type for the global grids.
"""
abstract type GridsT end



include("FEMGrid.jl")
include("SMGrid.jl")

export rfem_grid
export afem_grid
export SMGridDataT
export asm_grid
export periodic_grid
export mode_list
export ifft_size
export mode_label
export inst_grid



include("GlobalGrids.jl")


export GridsT
export ContGridsT
export FSSGridsT
export FFSGridsT
export FFFGridsT
export init_grids
export inst_grids



include("Problem.jl")

export ProblemT
export GeoParamsT
export FLRT
export init_problem


include("Output.jl")

export EvalsT
export find_ind


include("Mapping.jl")

export MapGridsT

end