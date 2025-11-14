#these are working now, save a little bit of allocations.
#could probably be combined into one thing!
#simple in place form , which is just passed on the apropriate shit being passed in.
#try generalise later!
function local_to_global!(x1::Array{Float64}, dx1::Array{Float64}, ind::Int64, ξ::Array{Float64}, grid::Array{Float64})

    dx1[1] = grid[ind+1] - grid[ind]

    #local coords are mapped from (-1, 1) to (0, 1)
    #then scaled by dx
    #finally, shiften to the proper place in the grid.
    @. x1 = dx1[1] * ((ξ+1)/2) + grid[ind]

    return dx1[1] / 2
end


function local_to_global!(x1::Array{Float64}, x2::Array{Float64}, dx::Array{Float64}, x1ind::Int64, x2ind::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen)

    dx[1] = x1grid[x1ind+1] - x1grid[x1ind]

    if x2ind == length(x2grid)
        dx[2] = 2π + x2grid[1] - x2grid[end]
    else
        dx[2] = x2grid[x2ind+1] - x2grid[x2ind]
    end

    #local coords are mapped from (-1, 1) to (0, 1)
    #then scaled by dx
    #finally, shiften to the proper place in the grid.
    @. x1 = dx[1] * ((ξx1+1)/2) + x1grid[x1ind]
    @. x2 = dx[2] * ((ξx2+1)/2) + x2grid[x2ind]

    return dx[1] * dx[2] / 4
end

#
function local_to_global!(x1::Array{Float64}, x2::Array{Float64}, x3::Array{Float64}, dx::Array{Float64}, x1ind::Int64, x2ind::Int64, x3ind::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, ξx3::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen, x3grid::StepRangeLen)

    dx[1] = x1grid[x1ind+1] - x1grid[x1ind]

    if x2ind == length(x2grid)
        dx[2] = 2π + x2grid[1] - x2grid[end]
    else
        dx[2] = x2grid[x2ind+1] - x2grid[x2ind]
    end

    if x3ind == length(x3grid)
        dx[3] = 2π + x3grid[1] - x3grid[end]
    else
        dx[3] = x3grid[x3ind+1] - x3grid[x3ind]
    end

    #local coords are mapped from (-1, 1) to (0, 1)
    #then scaled by dx
    #finally, shiften to the proper place in the grid.
    @. x1 = dx[1] * ((ξx1+1)/2) + x1grid[x1ind]
    @. x2 = dx[2] * ((ξx2+1)/2) + x2grid[x2ind]
    @. x3 = dx[3] * ((ξx3+1)/2) + x3grid[x3ind]

    return dx[1] * dx[2] * dx[3] / 8
end

