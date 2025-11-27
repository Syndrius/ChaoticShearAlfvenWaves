"""
    local_to_global!(x1::Array{Float64}, dx1::Array{Float64}, ind::Int64, ξ::Array{Float64}, grid::Array{Float64})

Converts the local grid used by finite elements to the appropriate values in the global grid.
"""
function local_to_global!(x1::Array{Float64}, dx1::Array{Float64}, ind::Int64, ξ::Array{Float64}, grid::Array{Float64})

    dx1[1] = grid[ind+1] - grid[ind]

    #local coords are mapped from (-1, 1) to (0, 1)
    #then scaled by dx
    #finally, shifted to the proper place in the grid.
    @. x1 = dx1[1] * ((ξ+1)/2) + grid[ind]

    return dx1[1] / 2
end


"""
    local_to_global!(x1::Array{Float64}, x2::Array{Float64}, dx::Array{Float64}, x1ind::Int64, x2ind::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen)

Converts the local grid used by finite elements to the appropriate values in the global grid.
"""
function local_to_global!(x1::Array{Float64}, x2::Array{Float64}, dx::Array{Float64}, x1ind::Int64, x2ind::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen)

    dx[1] = x1grid[x1ind+1] - x1grid[x1ind]

    if x2ind == length(x2grid)
        dx[2] = 2π + x2grid[1] - x2grid[end]
    else
        dx[2] = x2grid[x2ind+1] - x2grid[x2ind]
    end

    #local coords are mapped from (-1, 1) to (0, 1)
    #then scaled by dx
    #finally, shifted to the proper place in the grid.
    @. x1 = dx[1] * ((ξx1+1)/2) + x1grid[x1ind]
    @. x2 = dx[2] * ((ξx2+1)/2) + x2grid[x2ind]

    return dx[1] * dx[2] / 4
end

"""
    local_to_global!(x1::Array{Float64}, x2::Array{Float64}, x3::Array{Float64}, dx::Array{Float64}, x1ind::Int64, x2ind::Int64, x3ind::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, ξx3::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen, x3grid::StepRangeLen)

Converts the local grid used by finite elements to the appropriate values in the global grid.
"""

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
    #finally, shifted to the proper place in the grid.
    @. x1 = dx[1] * ((ξx1+1)/2) + x1grid[x1ind]
    @. x2 = dx[2] * ((ξx2+1)/2) + x2grid[x2ind]
    @. x3 = dx[3] * ((ξx3+1)/2) + x3grid[x3ind]

    return dx[1] * dx[2] * dx[3] / 8
end


"""
    global_to_local(ind::Int64, grid::AbstractArray{Float64}, x::Float64)

Maps the global grid into t alocal grid defined between -1 and 1.
Used for Interpolation.
"""
function global_to_local(ind::Int64, grid::AbstractArray{Float64}, x::Float64)
    #this is essentially an inverse of our local to global functions
    #this takes the global segment, x∈[x_i, x_{i+1}] to local domain, ξ∈[-1, 1]
    #guard in case point is exactly on grid
    if grid[ind]==x
        #extra edge cases only for derivatives.
        ξ = -1.0
        ind1 = ind
        #assumes 2π periodicity
        if ind1 == length(grid)
            ind2 = 1
            Δx = 2π + (grid[ind2] - grid[ind1])  #global difference
        else
            ind2 = ind+1 
            Δx = grid[ind2] - grid[ind] 
        end
        inds = [ind1, ind2]
    #periodic case
    elseif x > grid[end]
        inds = [length(grid), 1]
        Δx = 2π + (grid[1] - grid[end])
        ξ = (x - grid[end]) * 2 / Δx - 1
    elseif grid[ind] < x
        inds = [ind, ind+1]
        Δx = grid[ind+1] - grid[ind]
        ξ = (x - grid[ind]) * 2 / Δx - 1
    else
        #may need to add some other edge cases 
        #in particular, the qfm grid not being maximal can cause problems
        #typically this is easily fixed by making the mapped grid smaller.
        inds = [ind-1, ind]
        Δx = grid[ind] - grid[ind-1]
        ξ = (x - grid[ind-1]) * 2 / Δx - 1
    end
    return ξ, inds, Δx
end


