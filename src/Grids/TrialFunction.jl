
#stupid name, stupid location
#don't know if this is better or worse than in weakform.

function update_trial_function!(Φ::Array{ComplexF64, 3}, S::HB1dT, m::Int64, n::Int64, jac::Float64, ts::Array{Float64, 2})

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
    create_global_basis!(Φ::Array{Float64, 5}, S::HB2d, m::Int64, n::Int64, dx1::Float64, dx2::Float64, ts::Array{Float64, 5})

Converts the local basis functions into global, by mapping ξ∈[-1, 1] tp x∈[x_i, x_{i+1}].
Additionally, m represents a phase factor for setting the '0' of the fourier transform to target specific modes.
"""
function update_trial_function!(Φ::Array{ComplexF64, 5}, S::HB2dT, m::Int64, n::Int64, dx1::Float64, dx2::Float64, ts::Array{Float64, 4})

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
    create_global_basis!(Φ::Array{Float64, 7}, S::HB3d, m::Int64, n::Int64, dx1::Float64, dx2::Float64, dx3::Float64, ts::Array{Float64, 7})

Converts the local basis functions into global, by mapping ξ∈[-1, 1] tp x∈[x_i, x_{i+1}].
Additionally, m and n represent phase factors for setting the '0' of the fourier transform to target specific modes.
"""
function update_trial_function!(Φ::Array{ComplexF64, 7}, S::HB3dT, m::Int64, n::Int64, dx1::Float64, dx2::Float64, dx3::Float64, ts::Array{Float64, 6})


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


