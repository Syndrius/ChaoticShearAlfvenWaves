
"""
    interpolation(x1::Float64, x2::Float64, x3::Float64, ϕ::Array{ComplexF64, 4}, x1grid::Array{Float64}, x2grid::AbstractArray{Float64}, x3grid::AbstractArray{Float64})

Interpolates the eigenfunction Φ at (x1, x2, x3) based on the Hermite basis functions.
"""
function interpolation(x1::Float64, x2::Float64, x3::Float64, ϕ::Array{ComplexF64, 4}, x1grid::Array{Float64}, x2grid::AbstractArray{Float64}, x3grid::AbstractArray{Float64})
    #first find the index of the grid point closest to the target point
    x1ind = argmin(abs.(x1grid .- x1))
    x2ind = argmin(abs.(x2grid .- x2))
    x3ind = argmin(abs.(x3grid .- x3))

    #using this point, we find the two indices left and right of the point, the distance from the left point (ξ) and the disatnce between the two points.
    ξ1, inds1, Δx1 = global_to_local(x1ind, x1grid, x1)
    ξ2, inds2, Δx2 = global_to_local(x2ind, x2grid, x2)
    ξ3, inds3, Δx3 = global_to_local(x3ind, x3grid, x3)

    ϕ_int = 0.0+0.0im
    gid = grid_id
    bid = basis_id
    for h3 in 1:4, h2 in 1:4, h1 in 1:4
        gi = (inds1[gid[h1]+1], inds2[gid[h2]+1], inds3[gid[h3]+1])
        bi = 1 + 4*bid[h1] + 2*bid[h2] + 1*bid[h3] #cannot be sure this is actually valid yet!
        ϕ_int += ϕ[gi..., bi] * hb(ξ1, h1, Δx1) * hb(ξ2, h2, Δx2) * hb(ξ3, h3, Δx3)
    end
    return ϕ_int

end

"""
    interpolation(x1::Float64, x2::Float64, x3::Float64, ϕ::Array{ComplexF64, 4}, x1grid::Array{Float64}, x2grid::AbstractArray{Float64}, x3grid::AbstractArray{Float64}, x1d::Int64, x2d::Int64, x3d::Int64)

Interpolates the eigenfunction Φ for specified deriavtes in each dimension.
"""
function interpolation(x1::Float64, x2::Float64, x3::Float64, ϕ::Array{ComplexF64, 4}, x1grid::Array{Float64}, x2grid::AbstractArray{Float64}, x3grid::AbstractArray{Float64}, x1d::Int64, x2d::Int64, x3d::Int64)
    #first find the index of the grid point closest to the target point
    x1ind = argmin(abs.(x1grid .- x1))
    x2ind = argmin(abs.(x2grid .- x2))
    x3ind = argmin(abs.(x3grid .- x3))

    #using this point, we find the two indices left and right of the point, the distance from the left point (ξ) and the disatnce between the two points.
    ξ1, inds1, Δx1 = global_to_local(x1ind, x1grid, x1)
    ξ2, inds2, Δx2 = global_to_local(x2ind, x2grid, x2)
    ξ3, inds3, Δx3 = global_to_local(x3ind, x3grid, x3)

    ϕ_int = 0.0+0.0im
    gid = grid_id
    bid = basis_id
    for h1 in 1:4, h2 in 1:4, h3 in 1:4
        gi = (inds1[gid[h1]+1], inds2[gid[h2]+1], inds3[gid[h3]+1])
        bi = 1 + 4*bid[h1] + 2*bid[h2] + 1*bid[h3] #cannot be sure this is actually valid yet!
        ϕ_int += ϕ[gi..., bi] * hb(ξ1, h1, Δx1, x1d) * hb(ξ2, h2, Δx2, x2d) * hb(ξ3, h3, Δx3, x3d)
    end
    return ϕ_int

end


"""
    hb(t::Float64, h::Int64, dt::Float64, order::Int64=0)

Function that returns the appropriate hermite basis function based on h.
"""
function hb(ξ::Float64, h::Int64, Δx::Float64, order::Int64=0)
    #Δx is due to arbitrary interval, see wikipedia page.
    #change from ξ∈[-1, 1] to x∈[x_i, x_{i+1}]
    jac = Δx / 2 
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


