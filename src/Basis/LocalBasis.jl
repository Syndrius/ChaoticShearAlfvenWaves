
"""
    create_global_basis!(Φ::Array{ComplexF64, 3}, H::Array{Float64, 2}, dH::Array{Float64, 2}, ddH::Array{Float64, 2}, m::Int64, n::Int64, jac::Float64, ts::Array{Float64, 2})

Converts the local basis functions into global, by mapping ξ∈[-1, 1] tp x∈[x_i, x_{i+1}].

### Args
Φ::Array{ComplexF64, 3} - 9x4xgp matrix, stores the 9 derivtaives of the potential for each 4 Hermite basis functions at each gauss point.
H::Array{Float64, 2} - Hermite basis functions.
dH::Array{Float64, 2} - Derivtive of Hermite basis functions.
ddH::Array{Float64, 2} - Second derivtive of Hermite basis functions.
m::Int64 - Poloidal mode number, i.e. scale factor for derivative of Fourier basis.
n::Int64 - Toroidal mode number, i.e. scale factor for derivative of Fourier basis.
jac::Float64 - Jacobian of local to global transformation, not to be confused with geometric Jacobian.
ts::Array{Float64, 2} - Array that stores the tangent scale, relfecting the extra scaling the tangent basis functions, h10 and h11, need in the local to global transformation.
"""
function create_global_basis!(Φ::Array{ComplexF64, 3}, S::HB1d, m::Int64, n::Int64, jac::Float64, ts::Array{Float64, 2})

    #reset the tangent scaling to 1.
    ts .= 1.0

    #tangent basis functions need additional scaling when transforming from global ξ∈[-1, 1]
    #to global x∈[x_i, x_{i+1}]
    #to ensure they have the correct properties, i.e. derivative equals 1, at the boundaries.
    ts[2, :] *= jac
    ts[4, :] *= jac

    #Φ_x1
    @. Φ[:, 1, :] = S.dH / jac * ts
    #Φ_x2
    @. Φ[:, 2, :] = S.H * m * 1im * ts
    #Φ_x3
    @. Φ[:, 3, :] = S.H * n * 1im * ts
    #Φ_x1x1
    @. Φ[:, 4, :] = S.ddH / jac^2 * ts
    #Φ_x1x2
    @. Φ[:, 5, :] = S.dH * m * 1im / jac * ts
    #Φ_x1x3   
    @. Φ[:, 6, :] = S.dH * n * 1im / jac * ts
    #Φ_x2x2
    @. Φ[:, 7, :] = S.H * (-m^2) * ts
    #Φ_x2x3
    @. Φ[:, 8, :] = S.H * (-m*n) * ts
    #Φ_x3x3
    @. Φ[:, 9, :] = S.H * (-n^2) * ts

end


"""
    create_global_basis!(Φ::Array{ComplexF64, 5}, S::HB2d, m::Int64, n::Int64, dx1::Float64, dx2::Float64, ts::Array{Float64, 5})

Converts the local basis functions into global, by mapping ξ∈[-1, 1] tp x∈[x_i, x_{i+1}].
Additionally, m represents a phase factor for setting the '0' of the fourier transform to target specific modes.
"""
function create_global_basis!(Φ::Array{ComplexF64, 5}, S::HB2d, m::Int64, n::Int64, dx1::Float64, dx2::Float64, ts::Array{Float64, 4})

    #jacobian of the local to global transformation in each dimension
    x1jac = dx1 / 2
    x2jac = dx2 / 2

    #reset the tangent scaling to 1.
    ts .= 1.0

    #tangent basis functions need additional scaling when transforming from global ξ∈[-1, 1]
    #to global x∈[x_i, x_{i+1}]
    #to ensure they have the correct properties, i.e. derivative equals 1, at the boundaries.
    ts[2, :, :, :] *= x1jac
    ts[4, :, :, :] *= x1jac
    ts[:, 2, :, :] *= x2jac
    ts[:, 4, :, :] *= x2jac
    
    #Φ_x1
    @. Φ[:, :, 1, :, :] = S.dHx1 / x1jac * ts
    #Φ_x2
    @. Φ[:, :, 2, :, :] = (S.dHx2 / x2jac + S.H * 1im * m) * ts
    #Φ_x3
    @. Φ[:, :, 3, :, :] = S.H * n * 1im * ts
    #Φ_x1x1
    @. Φ[:, :, 4, :, :] = S.ddHx1x1 / x1jac^2 * ts
    #Φ_x1x2
    @. Φ[:, :, 5, :, :] = (S.ddHx1x2 / (x1jac * x2jac) + S.dHx1 * 1im * m / x1jac) * ts
    #Φ_x1x3   
    @. Φ[:, :, 6, :, :] = S.dHx1 * n * 1im / x1jac * ts
    #Φ_x2x2
    @. Φ[:, :, 7, :, :] = (S.ddHx2x2 / x2jac^2 + 2 * S.dHx2 * 1im * m / x2jac - m^2 * S.H) * ts
    #Φ_x2x3
    @. Φ[:, :, 8, :, :] = 1im * n * (S.dHx2 / x2jac + S.H * 1im * m) * ts
    #Φ_x3x3
    @. Φ[:, :, 9, :, :] = S.H * (-n^2) * ts

end




"""
    create_global_basis!(Φ::Array{ComplexF64, 7}, S::HB3d, m::Int64, n::Int64, dx1::Float64, dx2::Float64, dx3::Float64, ts::Array{Float64, 7})

Converts the local basis functions into global, by mapping ξ∈[-1, 1] tp x∈[x_i, x_{i+1}].
Additionally, m and n represent phase factors for setting the '0' of the fourier transform to target specific modes.
"""
function create_global_basis!(Φ::Array{ComplexF64, 7}, S::HB3d, m::Int64, n::Int64, dx1::Float64, dx2::Float64, dx3::Float64, ts::Array{Float64, 6})


    #jacobian of the local to global transformation in each dimension
    x1jac = dx1 / 2
    x2jac = dx2 / 2
    x3jac = dx3 / 2

    #reset the tangent scaling to 1.
    ts .= 1.0

    #tangent basis functions need additional scaling when transforming from global ξ∈[-1, 1]
    #to global x∈[x_i, x_{i+1}]
    #to ensure they have the correct properties, i.e. derivative equals 1, at the boundaries.
    ts[2, :, :, :, :, :] *= x1jac
    ts[4, :, :, :, :, :] *= x1jac
    ts[:, 2, :, :, :, :] *= x2jac
    ts[:, 4, :, :, :, :] *= x2jac
    ts[:, :, 2, :, :, :, :] *= x3jac
    ts[:, :, 4, :, :, :, :] *= x3jac
    
    #Φ_x1
    @. Φ[:, :, :, 1, :, :, :] = S.dHx1 / x1jac * ts
    #Φ_x2
    @. Φ[:, :, :, 2, :, :, :] = (S.dHx2 / x2jac + S.H * 1im * m) * ts
    #Φ_x3
    @. Φ[:, :, :, 3, :, :, :] = (S.dHx3 / x3jac + S.H * 1im * n) * ts
    #Φ_x1x1
    @. Φ[:, :, :, 4, :, :, :] = S.ddHx1x1 / x1jac^2 * ts
    #Φ_x1x2
    @. Φ[:, :, :, 5, :, :, :] = (S.ddHx1x2 / (x1jac * x2jac) + S.dHx1 * 1im * m / x1jac) * ts
    #Φ_x1x3   
    @. Φ[:, :, :, 6, :, :, :] = (S.ddHx1x3 / (x1jac * x3jac) + S.dHx1 * 1im * n / x1jac) * ts
    #Φ_x2x2
    @. Φ[:, :, :, 7, :, :, :] = (S.ddHx2x2 / x2jac^2 + 2 * S.dHx2 * 1im * m / x2jac - m^2 * S.H) * ts
    #Φ_x2x3
    @. Φ[:, :, :, 8, :, :, :] = (S.ddHx2x3 / (x2jac * x3jac) + S.dHx3 * 1im * m / x3jac + S.dHx2 * 1im * n / x2jac - m * n * S.H) * ts
    #Φ_x3x3
    @. Φ[:, :, :, 9, :, :, :] = (S.ddHx3x3 / x3jac^2 + 2 * S.dHx3 * 1im * n / x3jac - n^2 * S.H) * ts

end


"""
    local_to_global(node::Int64, ξ::Array{Float64}, grid::Array{Float64})

Converts the local grid where the finite element basis is defined to the global coordinates.
"""
function local_to_global(node::Int64, ξ::Array{Float64}, grid::Array{Float64})

    #TODO
    #this could be made non-allocating

    #should this be in basis??
    dx1 = grid[node+1] - grid[node]

    #first we map the (-1, 1) global coords to (0, 1)
    mp = @. (ξ+1) / 2

    #then we sale by the width of grid points
    mp = mp .* dx1

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    x1global = mp .+ grid[node]

    
    return x1global, dx1
end 




"""
    local_to_global(x1node::Int64, x2node::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen)

Converts the local grid where the finite element basis is defined to the global coordinates.
"""
function local_to_global(x1node::Int64, x2node::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen)

    #handles periodicity
    if x2node == length(x2grid)
        dx2 = 2π + x2grid[1] - x2grid[x2node]
    else
        dx2 = x2grid[x2node+1] - x2grid[x2node]
    end


    dx1 = x1grid[x1node+1] - x1grid[x1node]
    

    #first we map the (-1, 1) global coords to (0, 1)
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

Converts the local grid where the finite element basis is defined to the global coordinates.
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
    
    #first we map the (-1, 1) global coords to (0, 1)
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
