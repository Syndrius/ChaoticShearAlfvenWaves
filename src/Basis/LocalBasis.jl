
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

    #Φ_x1
    @. Φ[:, 1, :] = S.dH / jac
    #Φ_x2
    @. Φ[:, 2, :] = S.H * m * 1im
    #Φ_x3
    @. Φ[:, 3, :] = S.H * n * 1im
    #Φ_x1x1
    @. Φ[:, 4, :] = S.ddH / jac^2
    #Φ_x1x2
    @. Φ[:, 5, :] = S.dH * m * 1im / jac
    #Φ_x1x3   
    @. Φ[:, 6, :] = S.dH * n * 1im / jac
    #Φ_x2x2
    @. Φ[:, 7, :] = S.H * (-m^2)
    #Φ_x2x3
    @. Φ[:, 8, :] = S.H * (-m*n)
    #Φ_x3x3
    @. Φ[:, 9, :] = S.H * (-n^2)

end


"""
    create_local_basis!(Φ::Array{ComplexF64, 5}, S::Array{Float64, 4}, dSr::Array{Float64, 4}, dSx2::Array{Float64, 4}, ddSrr::Array{Float64, 4}, ddSrx2::Array{Float64, 4}, ddSx2x2::Array{Float64, 4}, m::Int64, n::Int64, dx1::Float64, dx2::Float64)

Modifies the matrix storing the basis functions, Φ, to reflect the local derivatives being taken. Includes a phase factor modification for x2.
"""
function create_local_basis!(Φ::Array{ComplexF64, 5}, S::HB2d, m::Int64, n::Int64, dx1::Float64, dx2::Float64)

    x1jac = dx1 / 2
    x2jac = dx2 / 2
    
    #Φ_x1
    @. Φ[:, :, 1, :, :] = S.dHx1 / x1jac 
    #Φ_x2
    @. Φ[:, :, 2, :, :] = S.dHx2 / x2jac + S.H * 1im * m
    #Φ_x3
    @. Φ[:, :, 3, :, :] = S.H * n * 1im 
    #Φ_x1x1
    @. Φ[:, :, 4, :, :] = S.ddHx1x1 / x1jac^2 
    #Φ_x1x2
    @. Φ[:, :, 5, :, :] = S.ddHx1x2 / (x1jac * x2jac) + S.dHx1 * 1im * m / x1jac
    #Φ_x1x3   
    @. Φ[:, :, 6, :, :] = S.dHx1 * n * 1im / x1jac 
    #Φ_x2x2
    @. Φ[:, :, 7, :, :] = S.ddHx2x2 / x2jac^2 + 2 * S.dHx2 * 1im * m / x2jac - m^2 * S.H
    #Φ_x2x3
    @. Φ[:, :, 8, :, :] = 1im * n * (S.dHx2 / x2jac + S.H * 1im * m)
    #Φ_x3x3
    @. Φ[:, :, 9, :, :] = S.H * (-n^2) 

end




"""
    create_local_basis!(Φ::Array{ComplexF64, 7}, S::Array{Float64, 6}, dSr::Array{Float64, 6}, dSx2::Array{Float64, 6}, dSx3::Array{Float64, 6}, ddSrr::Array{Float64, 6}, ddSrx2::Array{Float64, 6}, ddSrx3::Array{Float64, 6}, ddSx2x2::Array{Float64, 6}, ddSx2x3::Array{Float64, 6}, ddSx3x3::Array{Float64, 6}, m::Int64, n::Int64, dx1::Float64, dx2::Float64, dx3::Float64)

Modifies the matrix storing the basis functions, Φ, to reflect the local derivatives being taken. Includes a phase factor modification for x2 and x3.
"""
function create_local_basis!(Φ::Array{ComplexF64, 7}, S::HB3d, m::Int64, n::Int64, dx1::Float64, dx2::Float64, dx3::Float64)


    x1jac = dx1 / 2
    x2jac = dx2 / 2
    x3jac = dx3 / 2
    
    #Φ_x1
    @. Φ[:, :, :, 1, :, :, :] = S.dHx1 / x1jac 
    #Φ_x2
    @. Φ[:, :, :, 2, :, :, :] = S.dHx2 / x2jac + S.H * 1im * m
    #Φ_x3
    @. Φ[:, :, :, 3, :, :, :] = S.dHx3 / x3jac + S.H * 1im * n
    #Φ_x1x1
    @. Φ[:, :, :, 4, :, :, :] = S.ddHx1x1 / x1jac^2 
    #Φ_x1x2
    @. Φ[:, :, :, 5, :, :, :] = S.ddHx1x2 / (x1jac * x2jac) + S.dHx1 * 1im * m / x1jac
    #Φ_x1x3   
    @. Φ[:, :, :, 6, :, :, :] = S.ddHx1x3 / (x1jac * x3jac) + S.dHx1 * 1im * n / x1jac 
    #Φ_x2x2
    @. Φ[:, :, :, 7, :, :, :] = S.ddHx2x2 / x2jac^2 + 2 * S.dHx2 * 1im * m / x2jac - m^2 * S.H
    #Φ_x2x3
    @. Φ[:, :, :, 8, :, :, :] = S.ddHx2x3 / (x2jac * x3jac) + S.dHx3 * 1im * m / x3jac + S.dHx2 * 1im * n / x2jac - m * n * S.H
    #Φ_x3x3
    @. Φ[:, :, :, 9, :, :, :] = S.ddHx3x3 / x3jac^2 + 2 * S.dHx3 * 1im * n / x3jac - n^2 * S.H

end


"""
    local_to_global(node::Int64, ξ::Array{Float64}, grid::Array{Float64})

Converts the local grid where the finite elements are defined to the global coordinates.
"""
function local_to_global(node::Int64, ξ::Array{Float64}, grid::Array{Float64})

    #should this be in basis??
    dx1 = grid[node+1] - grid[node]

    #first we map the (-1, 1) local coords to (0, 1)
    mp = @. (ξ+1) / 2

    #then we sale by the width of grid points
    mp = mp .* dx1

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    x1global = mp .+ grid[node]

    
    return x1global, dx1
end 




"""
    function local_to_global(x1node::Int64, x2node::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen)

Converts the local grid where the finite elements are defined to the global coordinates.
"""
function local_to_global(x1node::Int64, x2node::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen)

    #handles periodicity
    if x2node == length(x2grid)
        dx2 = 2π + x2grid[1] - x2grid[x2node]
    else
        dx2 = x2grid[x2node+1] - x2grid[x2node]
    end


    dx1 = x1grid[x1node+1] - x1grid[x1node]
    

    #first we map the (-1, 1) local coords to (0, 1)
    mpx1 = @. (ξx1+1) / 2
    mpx2 = @. (ξx2+1) / 2

    #then we sale by the width of grid points
    mpx1 = mpx1 .* dx1
    mpx2 = mpx2 .* dx2

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    x1global = mpx1 .+ x1grid[x1node]
    x2global = mpx2 .+ x2grid[x2node]

    
    return x1global, x2global, dx1, dx2
end 


"""
    local_to_global(x1node::Int64, x2node::Int64, x3node::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, ξx3::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen, x3grid::StepRangeLen)

Converts the local grid where the finite elements are defined to the global coordinates.
"""
function local_to_global(x1node::Int64, x2node::Int64, x3node::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, ξx3::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen, x3grid::StepRangeLen)

    #handles periodicity
    if x2node == length(x2grid)
        dx2 = 2π + x2grid[1] - x2grid[x2node]
    else
        dx2 = x2grid[x2node+1] - x2grid[x2node]
    end

    if x3node == length(x3grid)
        dx3 = 2π + x3grid[1] - x3grid[x3node]
    else
        dx3 = x3grid[x3node+1] - x3grid[x3node]
    end

    dx1 = x1grid[x1node+1] - x1grid[x1node]
    
    #first we map the (-1, 1) local coords to (0, 1)
    mpx1 = @. (ξx1+1) / 2
    mpx2 = @. (ξx2+1) / 2
    mpx3 = @. (ξx3+1) / 2

    #then we sale by the width of grid points
    mpx1 = mpx1 .* dx1
    mpx2 = mpx2 .* dx2
    mpx3 = mpx3 .* dx3

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    x1global = mpx1 .+ x1grid[x1node]
    x2global = mpx2 .+ x2grid[x2node]
    x3global = mpx3 .+ x3grid[x3node]

    
    return x1global, x2global, x3global, dx1, dx2, dx3
end 
