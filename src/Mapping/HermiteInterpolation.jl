

#currently using in build interpolation instead.
#this should be more accurate though!
function hermite_interpolation(x1::Float64, x2::Float64, x3::Float64, ϕ::Array{ComplexF64, 4}, x1grid::Array{Float64}, x2grid::AbstractArray{Float64}, x3grid::AbstractArray{Float64})
    #first find the index of the grid point closest to the target point
    x1ind = find_ind(x1grid, x1)
    x2ind = find_ind(x2grid, x2)
    x3ind = find_ind(x3grid, x3)

    #using theis point, we find the two indices left and right of the point, the distance from the left point (ξ) and the disatnce between the two points.
    ξ1, inds1, dx1 = global_to_local(x1ind, x1grid, x1)
    ξ2, inds2, dx2 = global_to_local(x2ind, x2grid, x2)
    ξ3, inds3, dx3 = global_to_local(x3ind, x3grid, x3)

    #going to assume the grid_id/basis_id is already defined. May need a more solid import from basis or whatever.

    ϕ_int = 0.0+0.0im
    gid = Indexing.grid_id
    bid = Indexing.basis_id
    for h1 in 1:4, h2 in 1:4, h3 in 1:4
        gi = (inds1[gid[h1]+1], inds2[gid[h2]+1], inds3[gid[h3]+1])
        bi = 1 + bid[h1] + 2*bid[h2] + 4*bid[h3] #cannot be sure this is actually valid yet!
        ϕ_int += ϕ[gi..., bi] * hb(ξ1, h1, dx1) * hb(ξ2, h2, dx2) * hb(ξ3, h3, dx3)
    end
    return ϕ_int

end

#forces periodicity, assumes the grids have the usual structure.
#for r/κ chosen interpoaltion grid should never go beyond real grid. If so, it will be treated as periodic.
#may cause a problemo for mapping island to outside.
#but this should just be set to zero elsewhere.
function global_to_local(ind::Int64, grid::AbstractArray{Float64}, x::Float64)

    #guard in case point is exactly on grid
    if grid[ind]==x
        ξ = 0.0
        inds = [ind, ind]
        dx = 1.0 #choosing 1.0 due to division later, for this case this value should not matter. as Shape functions should evaluate to zero.
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
    #return (1+2*t)*(1-t)^2
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

function hb(t::Float64, h::Int64, dt::Float64)
    #additional jacobian term is very awkward.

    #dt is required as we are transforming from [-1, 1] to [0, 1]
    #and from global to local, so the derivatives need an extra jacobian.
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
