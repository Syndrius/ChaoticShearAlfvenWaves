"""
Struct storing the data required for computing the continuum. Simplified RFEMGridT, used for multiple dispatch.

### Fields
- N::Int64 - Number of grid points.
- start :: Real - Starting point of the grid.
- stop :: Real - Stopping point of the grid.
"""
@kwdef struct ContGridDataT <: GridDataT
    N :: Int64
    start :: Real = 0 
    stop :: Real = 1 
end


"""
    inst_grid(grid::ContGridDataT)

Instantiates the grid used for computation.
"""
function inst_grid(grid::ContGridDataT)
    #to prevent r=0 in the continuum case
    if grid.start == 0
        return LinRange(grid.start, grid.stop, grid.N+1)[2:end]
    else 
        return LinRange(grid.start, grid.stop, grid.N)
    end
end
