
"""
    create_local_basis!(Φ::Array{ComplexF64, 3}, H::Array{Float64, 2}, dH::Array{Float64, 2}, ddH::Array{Float64, 2}, m::Int64, n::Int64, jac::Float64)

Modifies the matrix storing the basis functions, Φ, to reflect the local derivatives being taken.

### Args
Φ::Array{ComplexF64, 3} - 9x4xgp matrix, stores the 9 derivtaives of the potential for each 4 Hermite basis functions at each gauss point.
H::Array{Float64, 2} - Hermite basis functions.
dH::Array{Float64, 2} - Derivtive of Hermite basis functions.
ddH::Array{Float64, 2} - Second derivtive of Hermite basis functions.
m::Int64 - Poloidal mode number, i.e. scale factor for derivative of Fourier basis.
n::Int64 - Toroidal mode number, i.e. scale factor for derivative of Fourier basis.
jac::Float64 - Jacobian of local to global transformation, not to be confused with geometric Jacobian.
"""
function create_local_basis!(Φ::Array{ComplexF64, 3}, S::HB1d, m::Int64, n::Int64, jac::Float64)

    #the middle index is chosen to increase efficiency with numerical integration.

    #Φ_s
    @. Φ[:, 1, :] = S.dH / jac
    #Φ_θ
    @. Φ[:, 2, :] = S.H * m * 1im
    #Φ_ζ
    @. Φ[:, 3, :] = S.H * n * 1im
    #Φ_ss
    @. Φ[:, 4, :] = S.ddH / jac^2
    #Φ_sθ
    @. Φ[:, 5, :] = S.dH * m * 1im / jac
    #Φ_sζ   
    @. Φ[:, 6, :] = S.dH * n * 1im / jac
    #Φ_θθ
    @. Φ[:, 7, :] = S.H * (-m^2)
    #Φ_θζ
    @. Φ[:, 8, :] = S.H * (-m*n)
    #Φ_ζζ
    @. Φ[:, 9, :] = S.H * (-n^2)

end


"""
    create_local_basis!(Φ::Array{ComplexF64, 5}, S::Array{Float64, 4}, dSr::Array{Float64, 4}, dSθ::Array{Float64, 4}, ddSrr::Array{Float64, 4}, ddSrθ::Array{Float64, 4}, ddSθθ::Array{Float64, 4}, m::Int64, n::Int64, dr::Float64, dθ::Float64)

Modifies the matrix storing the basis functions, Φ, to reflect the local derivatives being taken. Includes a phase factor modification for θ.
"""
function create_local_basis!(Φ::Array{ComplexF64, 5}, S::HB2d, m::Int64, n::Int64, dr::Float64, dθ::Float64)

    rjac = dr / 2
    θjac = dθ / 2
    
    #Φ_s
    @. Φ[:, :, 1, :, :] = S.dHr / rjac 
    #Φ_θ
    @. Φ[:, :, 2, :, :] = S.dHθ / θjac + S.H * 1im * m
    #Φ_ζ
    @. Φ[:, :, 3, :, :] = S.H * n * 1im 
    #Φ_ss
    @. Φ[:, :, 4, :, :] = S.ddHrr / rjac^2 
    #Φ_sθ
    @. Φ[:, :, 5, :, :] = S.ddHrθ / (rjac * θjac) + S.dHr * 1im * m / rjac
    #Φ_sζ   
    @. Φ[:, :, 6, :, :] = S.dHr * n * 1im / rjac 
    #Φ_θθ
    @. Φ[:, :, 7, :, :] = S.ddHθθ / θjac^2 + 2 * S.dHθ * 1im * m / θjac - m^2 * S.H
    #Φ_θζ
    @. Φ[:, :, 8, :, :] = 1im * n * (S.dHθ / θjac + S.H * 1im * m)
    #Φ_ζζ
    @. Φ[:, :, 9, :, :] = S.H * (-n^2) 

end




"""
    create_local_basis!(Φ::Array{ComplexF64, 7}, S::Array{Float64, 6}, dSr::Array{Float64, 6}, dSθ::Array{Float64, 6}, dSζ::Array{Float64, 6}, ddSrr::Array{Float64, 6}, ddSrθ::Array{Float64, 6}, ddSrζ::Array{Float64, 6}, ddSθθ::Array{Float64, 6}, ddSθζ::Array{Float64, 6}, ddSζζ::Array{Float64, 6}, m::Int64, n::Int64, dr::Float64, dθ::Float64, dζ::Float64)

Modifies the matrix storing the basis functions, Φ, to reflect the local derivatives being taken. Includes a phase factor modification for θ and ζ.
"""
function create_local_basis!(Φ::Array{ComplexF64, 7}, S::HB3d, m::Int64, n::Int64, dr::Float64, dθ::Float64, dζ::Float64)


    rjac = dr / 2
    θjac = dθ / 2
    ζjac = dζ / 2
    
    #Φ_s
    @. Φ[:, :, :, 1, :, :, :] = S.dHr / rjac 
    #Φ_θ
    @. Φ[:, :, :, 2, :, :, :] = S.dHθ / θjac + S.H * 1im * m
    #Φ_ζ
    @. Φ[:, :, :, 3, :, :, :] = S.dHζ / ζjac + S.H * 1im * n
    #Φ_ss
    @. Φ[:, :, :, 4, :, :, :] = S.ddHrr / rjac^2 
    #Φ_sθ
    @. Φ[:, :, :, 5, :, :, :] = S.ddHrθ / (rjac * θjac) + S.dHr * 1im * m / rjac
    #Φ_sζ   
    @. Φ[:, :, :, 6, :, :, :] = S.ddHrζ / (rjac * ζjac) + S.dHr * 1im * n / rjac 
    #Φ_θθ
    @. Φ[:, :, :, 7, :, :, :] = S.ddHθθ / θjac^2 + 2 * S.dHθ * 1im * m / θjac - m^2 * S.H
    #Φ_θζ
    @. Φ[:, :, :, 8, :, :, :] = S.ddHθζ / (θjac * ζjac) + S.dHζ * 1im * m / ζjac + S.dHθ * 1im * n / θjac - m * n * S.H
    #Φ_ζζ
    @. Φ[:, :, :, 9, :, :, :] = S.ddHζζ / ζjac^2 + 2 * S.dHζ * 1im * n / ζjac - n^2 * S.H

end


"""
    local_to_global(node::Int64, ξ::Array{Float64}, grid::Array{Float64})

Converts the local grid where the finite elements are defined to the global coordinates.
"""
function local_to_global(node::Int64, ξ::Array{Float64}, grid::Array{Float64})

    #should this be in basis??
    dr = grid[node+1] - grid[node]

    #first we map the (-1, 1) local coords to (0, 1)
    mp = @. (ξ+1) / 2

    #then we sale by the width of grid points
    mp = mp .* dr

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    rglobal = mp .+ grid[node]

    
    return rglobal, dr
end 




"""
    function local_to_global(rnode::Int64, θnode::Int64, ξr::Array{Float64}, ξθ::Array{Float64}, rgrid::Array{Float64}, θgrid::StepRangeLen)

Converts the local grid where the finite elements are defined to the global coordinates.
"""
function local_to_global(rnode::Int64, θnode::Int64, ξr::Array{Float64}, ξθ::Array{Float64}, rgrid::Array{Float64}, θgrid::StepRangeLen)

    #handles periodicity
    if θnode == length(θgrid)
        dθ = 2π + θgrid[1] - θgrid[θnode]
    else
        dθ = θgrid[θnode+1] - θgrid[θnode]
    end


    dr = rgrid[rnode+1] - rgrid[rnode]
    

    #first we map the (-1, 1) local coords to (0, 1)
    mpr = @. (ξr+1) / 2
    mpθ = @. (ξθ+1) / 2

    #then we sale by the width of grid points
    mpr = mpr .* dr
    mpθ = mpθ .* dθ

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    rglobal = mpr .+ rgrid[rnode]
    θglobal = mpθ .+ θgrid[θnode]

    
    return rglobal, θglobal, dr, dθ
end 


"""
    local_to_global(rnode::Int64, θnode::Int64, ζnode::Int64, ξr::Array{Float64}, ξθ::Array{Float64}, ξζ::Array{Float64}, rgrid::Array{Float64}, θgrid::StepRangeLen, ζgrid::StepRangeLen)

Converts the local grid where the finite elements are defined to the global coordinates.
"""
function local_to_global(rnode::Int64, θnode::Int64, ζnode::Int64, ξr::Array{Float64}, ξθ::Array{Float64}, ξζ::Array{Float64}, rgrid::Array{Float64}, θgrid::StepRangeLen, ζgrid::StepRangeLen)

    #handles periodicity
    if θnode == length(θgrid)
        dθ = 2π + θgrid[1] - θgrid[θnode]
    else
        dθ = θgrid[θnode+1] - θgrid[θnode]
    end

    if ζnode == length(ζgrid)
        dζ = 2π + ζgrid[1] - ζgrid[ζnode]
    else
        dζ = ζgrid[ζnode+1] - ζgrid[ζnode]
    end

    dr = rgrid[rnode+1] - rgrid[rnode]
    
    #first we map the (-1, 1) local coords to (0, 1)
    mpr = @. (ξr+1) / 2
    mpθ = @. (ξθ+1) / 2
    mpζ = @. (ξζ+1) / 2

    #then we sale by the width of grid points
    mpr = mpr .* dr
    mpθ = mpθ .* dθ
    mpζ = mpζ .* dζ

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    rglobal = mpr .+ rgrid[rnode]
    θglobal = mpθ .+ θgrid[θnode]
    ζglobal = mpζ .+ ζgrid[ζnode]

    
    return rglobal, θglobal, ζglobal, dr, dθ, dζ
end 
