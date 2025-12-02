"""
Struct storing the data required for computing the continuum. Simplified RadialFEMGridT, used for multiple dispatch.

### Fields
- N::Int64 - Number of grid points.
- start :: Real - Starting point of the grid.
- stop :: Real - Stopping point of the grid.
"""
struct ContGridT <: GridT
    N :: Int64
    start :: Float64
    stop :: Float64
end


"""
    inst_grid(grid::ContGridT)

Instantiates the grid used for computation.
"""
function inst_grid(grid::ContGridT)
    #to prevent r=0 in the continuum case
    if grid.start == 0
        return LinRange(grid.start, grid.stop, grid.N+1)[2:end]
    else 
        return LinRange(grid.start, grid.stop, grid.N)
    end
end
