#TODO - unsure if we ever actually need this!
#also the basis functions are perhaps wrong now
#and should be called from Basis.

#also does the derivatives 
#unsure if this really needs to be a different function.
function hermite_interpolation(x1::Float64, x2::Float64, x3::Float64, ϕ::Array{ComplexF64, 4}, x1grid::Array{Float64}, x2grid::AbstractArray{Float64}, x3grid::AbstractArray{Float64}, x1d::Int64, x2d::Int64, x3d::Int64)
    #first find the index of the grid point closest to the target point
    x1ind = find_ind(x1grid, x1)
    x2ind = find_ind(x2grid, x2)
    x3ind = find_ind(x3grid, x3)

    #using this point, we find the two indices left and right of the point, the distance from the left point (ξ) and the disatnce between the two points.
    ξ1, inds1, dx1 = global_to_local(x1ind, x1grid, x1)
    ξ2, inds2, dx2 = global_to_local(x2ind, x2grid, x2)
    ξ3, inds3, dx3 = global_to_local(x3ind, x3grid, x3)


    ϕ_int = 0.0+0.0im
    gid = Indexing.grid_id
    bid = Indexing.basis_id
    for h1 in 1:4, h2 in 1:4, h3 in 1:4
        gi = (inds1[gid[h1]+1], inds2[gid[h2]+1], inds3[gid[h3]+1])
        bi = 1 + 4*bid[h1] + 2*bid[h2] + 1*bid[h3] #cannot be sure this is actually valid yet!
        ϕ_int += ϕ[gi..., bi] * hb(ξ1, h1, dx1, x1d) * hb(ξ2, h2, dx2, x2d) * hb(ξ3, h3, dx3, x3d)
    end
    return ϕ_int

end


#unsure if this will actually work!
function hermite_interpolation(x1::Float64, ϕ::Array{ComplexF64, 2}, x1grid::Array{Float64}, x1d::Int64)
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
        ϕ_int += ϕ[gi, bi] * hb(ξ1, h1, dx1, x1d) 
    end
    return ϕ_int

end

"""
    hb(t::Float64, h::Int64, dt::Float64)

Function that returns the appropriate hermite basis function based on h.
"""
function hb(ξ::Float64, h::Int64, Δx::Float64, order::Int64)
    #Δx is due to arbitrary interval, see wikipedia page.
    #change from ξ∈[-1, 1] to x∈[x_i, x_{i+1}]
    jac = Δx / 2 #may need to be upside down!
    if order == 0
        if h==1
            return Basis.h00(ξ)
        elseif h==2
            return Basis.h10(ξ) * jac
        elseif h==3
            return Basis.h01(ξ)
        else
            return Basis.h11(ξ) * jac
        end
    elseif order == 1
        if h==1
            return Basis.dh00(ξ) / jac
        elseif h==2
            return Basis.dh10(ξ) 
        elseif h==3
            return Basis.dh01(ξ) / jac
        else
            return Basis.dh11(ξ) 
        end
    elseif order == 2
        if h==1
            return Basis.ddh00(ξ) / jac^2
        elseif h==2
            return Basis.ddh10(ξ) / jac
        elseif h==3
            return Basis.ddh01(ξ) / jac^2
        else
            return Basis.ddh11(ξ) / jac
        end
    end
end


"""
    hb(t::Float64, h::Int64, dt::Float64)

Function that returns the appropriate hermite basis function based on h.
"""
#=
function hb(t::Float64, h::Int64, dt::Float64, deriv::Int64)
    #additional jacobian term is very awkward.

    if h==1
        if deriv == 1
            return dh00(t) / dt
        elseif deriv == 2
            return ddh00(t) / dt^2
        else
            return h00(t)
        end
    elseif h==2
        if deriv == 1
            return dh10(t) 
        elseif deriv == 2
            return ddh10(t) / dt
        else
            return h10(t) * dt
        end
    elseif h==3
        if deriv == 1
            return dh01(t) / dt
        elseif deriv == 2
            return ddh01(t) / dt^2
        else
            return h01(t)
        end
    else
        if deriv == 1
            return dh11(t) 
        elseif deriv == 2
            return ddh11(t) / dt
        else
            return h11(t) * dt
        end
    end
end
=#
#derivatives of the hermite basis functions.
#=
function dh00(t)

    return 6*t^2 - 6*t
end

function dh10(t)
    return (1-t)^2 - 2*t*(1-t)
end

function dh01(t)
    return -6*t^2 + 6*t
end

function dh11(t)
    return 2*t*(t-1) + t^2
end

function ddh00(t)
    return 12*t - 6
end

function ddh10(t)
    return -4*(1-t)+2*t
end

function ddh01(t)
    return -12*t + 6
end

function ddh11(t)
    return 4*t + 2*(t-1)
end
=#
