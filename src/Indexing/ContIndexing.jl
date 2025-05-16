"""
    grid_to_index(x2ind::Int64, x3ind::Int64, grids::ContGridsT)

Converts the index of the 2d grid into the proper place in the matrix.
"""
function grid_to_index(x2ind::Int64, x3ind::Int64, grids::ContGridsT)

    Nn = grids.x3.N

    return 1 + (x3ind - 1) + (x2ind-1) * Nn

end


#this doesn't really work because evals are returned in order not based on the matrix structure
#this only really works for evals
function index_to_grid(i::Int64, grids::ContGridsT)

    #guess!
    x2 = mod(div(i-1, grids.x3.N), grids.x2.N) + 1
    x3 = mod(i-1, grids.x3.N) + 1

    return x2, x3

end

