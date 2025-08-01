"""
    hermite_interpolation(x1::Float64, x2::Float64, x3::Float64, ϕ::Array{ComplexF64, 4}, x1grid::Array{Float64}, x2grid::AbstractArray{Float64}, x3grid::AbstractArray{Float64})

Interpolates the eigenfunction ϕ at (x1, x2, x3) based on the Hermite basis functions.
"""
function hermite_interpolation(x1::Float64, x2::Float64, x3::Float64, ϕ::Array{ComplexF64, 4}, x1grid::Array{Float64}, x2grid::AbstractArray{Float64}, x3grid::AbstractArray{Float64})
    #first find the index of the grid point closest to the target point
    x1ind = find_ind(x1grid, x1)
    x2ind = find_ind(x2grid, x2)
    x3ind = find_ind(x3grid, x3)

    #using this point, we find the two indices left and right of the point, the distance from the left point (ξ) and the disatnce between the two points.
    ξ1, inds1, dx1 = global_to_local(x1ind, x1grid, x1)
    ξ2, inds2, dx2 = global_to_local(x2ind, x2grid, x2)
    ξ3, inds3, dx3 = global_to_local(x3ind, x3grid, x3)


    #so something is wrong here.
    ϕ_int = 0.0+0.0im
    gid = Indexing.grid_id
    bid = Indexing.basis_id
    for h3 in 1:4, h2 in 1:4, h1 in 1:4
        gi = (inds1[gid[h1]+1], inds2[gid[h2]+1], inds3[gid[h3]+1])
        bi = 1 + 4*bid[h1] + 2*bid[h2] + 1*bid[h3] #cannot be sure this is actually valid yet!
        ϕ_int += ϕ[gi..., bi] * hb(ξ1, h1, dx1) * hb(ξ2, h2, dx2) * hb(ξ3, h3, dx3)
    end
    return ϕ_int

end


function hermite_interpolation(x1::Float64, ϕ::Array{ComplexF64, 2}, x1grid::Array{Float64})
    #first find the index of the grid point closest to the target point
    x1ind = find_ind(x1grid, x1)

    #using this point, we find the two indices left and right of the point, the distance from the left point (ξ) and the disatnce between the two points.
    ξ1, inds1, dx1 = global_to_local(x1ind, x1grid, x1)


    #so something is wrong here.
    ϕ_int = 0.0+0.0im
    gid = Indexing.grid_id
    bid = Indexing.basis_id
    for h1 in 1:4
        gi = inds1[gid[h1]+1]
        bi = 1 + bid[h1] 
        ϕ_int += ϕ[gi, bi] * hb(ξ1, h1, dx1) 
    end
    return ϕ_int

end


"""
    global_to_local(ind::Int64, grid::AbstractArray{Float64}, x::Float64)

Maps the global grid into t alocal grid defined between 0 and 1.
"""
function global_to_local(ind::Int64, grid::AbstractArray{Float64}, x::Float64)

    #guard in case point is exactly on grid
    if grid[ind]==x
        #extra edge cases only for derivatives.
        ξ = 0.0
        ind1 = ind
        #assumes 2π periodicity!
        if ind1 == length(grid)
            ind2 = 1
            dx = 2π + (grid[ind2] - grid[ind1])  #global difference
        else
            ind2 = ind+1 #doesn't matter
            dx = grid[ind2] - grid[ind] 
        end
        inds = [ind1, ind2]
    #periodic case
    elseif x > grid[end]
        inds = [length(grid), 1]
        dx = 2π + (grid[1] - grid[end])
        ξ = (x - grid[end]) / dx
    elseif grid[ind] < x
        inds = [ind, ind+1]
        dx = grid[ind+1] - grid[ind]
        ξ = (x - grid[ind]) / dx
    else
        #may need to add some other edge cases 
        #in particular, the qfm grid not being maximal can cause problemos
        #typically this is easily fixed by making the mapped grid smaller.
        inds = [ind-1, ind]
        dx = grid[ind] - grid[ind-1]
        ξ = (x - grid[ind-1]) / dx
    end
    return ξ, inds, dx
end


#needs to be moved
function h00(t::Float64)

    return 2*t^3 - 3*t^2 + 1
end

function h10(t::Float64)
    #return 2 * (t^3-2*t^2+t)
    return t*(1-t)^2
end

function h01(t::Float64)
    return -2t^3+3t^2
    #return t^2*(3-2*t)
end

function h11(t::Float64)
    #return 2*(t^3-t^2)
    return t^2*(t-1)
end

"""
    hb(t::Float64, h::Int64, dt::Float64)

Function that returns the appropriate hermite basis function based on h.
"""
function hb(t::Float64, h::Int64, dt::Float64)
    #additional jacobian term is very awkward.

    #dt is due to arbitrary interval, see wikipedia page.
    if h==1
        return h00(t)
    elseif h==2
        return h10(t) * dt
    elseif h==3
        return h01(t)
    else
        return h11(t) * dt
    end
end
