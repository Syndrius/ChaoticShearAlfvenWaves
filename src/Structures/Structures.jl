"""

Contains the structures and constructors used, includes
 - The ProblemT type.
 - The GridsT type.
 - The SolverT type.
 - The EvalsT type.

"""
module Structures


using ..Geometry
using ..Equilibrium

using FunctionWrappers
import FunctionWrappers: FunctionWrapper


"""
Abstract type for each grid.
"""
abstract type GridDataT end


"""
Abstract type for the global grids.
"""
abstract type GridsT end


"""
Abstract type for the Problem.
"""
abstract type ProblemT end


"""
Abstract type for Solving.
"""
abstract type SolverT end


include("FEMGrid.jl")
include("SMGrid.jl")
include("ContGrid.jl")

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
export init_grid
export init_geo
export inst_grids



include("FLR.jl")

export FLRT



include("Geometry.jl")

export GeoParamsT



include("Problem.jl")

export ProblemT
export init_problem
export IslProblemT
export TorProblemT #this is a bad name for this, given this can work with any metric....



include("Solver.jl")

export SolverT
export FullSpectrumSolverT
export ShiftInvertSolverT
export SliceSolverT
export init_solver

    

include("Output.jl")

export EvalsT
export find_ind

end
