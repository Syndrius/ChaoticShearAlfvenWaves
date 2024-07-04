

function grid_to_index(θind::Int64, ζind::Int64, grids::ContGridsT)

    Nn = grids.ζ.count

    return 1 + (ζind - 1) + (θind-1) * Nn

end


#this doesn't really work because evals are returned in order not based on the matrix structure
#this only really works for evals
function index_to_grid(i::Int64, grids::ContGridsT)

    #guess!
    θ = mod(div(i-1, grids.ζ.count), grids.θ.count) + 1
    ζ = mod(i-1, grids.ζ.count) + 1

    return θ, ζ

end


#still used!
function matrix_dim(grids::ContGridsT)

    return grids.θ.count * grids.ζ.count
end
