"""
    grid_to_index(θind::Int64, ζind::Int64, grids::ContGridsT)

Converts the index of the 2d grid into the proper place in the matrix.
"""
function grid_to_index(θind::Int64, ζind::Int64, grids::ContGridsT)

    Nn = grids.ζ.N

    return 1 + (ζind - 1) + (θind-1) * Nn

end


#this doesn't really work because evals are returned in order not based on the matrix structure
#this only really works for evals
function index_to_grid(i::Int64, grids::ContGridsT)

    #guess!
    θ = mod(div(i-1, grids.ζ.N), grids.θ.N) + 1
    ζ = mod(i-1, grids.ζ.N) + 1

    return θ, ζ

end

