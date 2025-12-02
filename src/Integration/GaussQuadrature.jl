
"""
    gauss_points(grids::FSSGridsT)

Returns the Gaussian integration points and corresponding weights
"""
function gauss_points(grids::FSSGridsT)

    return gausslegendre(grids.x1.gp)
end

"""
    gauss_points(grids::FFSGridsT)

Returns the Gaussian integration points and corresponding weights
"""
function gauss_points(grids::FFSGridsT)
    ξ1, wg1 = gausslegendre(grids.x1.gp)
    ξ2, wg2 = gausslegendre(grids.x2.gp)

    return ξ1, ξ2, wg1, wg2
end

"""
    gauss_points(grids::FFFGridsT)

Returns the Gaussian integration points and corresponding weights
"""
function gauss_points(grids::FFFGridsT)
    ξ1, wg1 = gausslegendre(grids.x1.gp)
    ξ2, wg2 = gausslegendre(grids.x2.gp)
    ξ3, wg3 = gausslegendre(grids.x3.gp)

    return ξ1, ξ2, ξ3, wg1, wg2, wg3
end


"""
    gauss_integrate(Ψ::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, Φ::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, mat::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, wg::Array{Float64}, jac::Float64, ngp::Int64)

Computes the numerical integration required for finite elements via ∫f(x)dx ≈ ∑w_i f(x_i) in 1d.
"""
function gauss_integrate(Ψ::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, Φ::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, mat::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, wg::Array{Float64}, jac::Float64, ngp::Int64)

    res = 0.0 + 0.0im
    for k in 1:ngp
    
        scale = wg[k] * jac
        for j in 1:9, i in 1:9
            
            res += @inbounds Ψ[i, k] * Φ[j, k] * mat[i, j, k] * scale

        end
    end
    return res 
end



"""
    gauss_integrate(Ψ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, Φ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, mat::SubArray{ComplexF64, 4, Array{ComplexF64, 5}}, wgx1::Array{Float64}, wgx2::Array{Float64}, jac::Float64, x1gp::Int64, x2gp::Int64)
    
Computes the numerical integration required for finite elements via ∫f(x)dx ≈ ∑w_i f(x_i) in 2d.
"""
function gauss_integrate(Ψ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, Φ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, mat::SubArray{ComplexF64, 4, Array{ComplexF64, 5}}, wgx1::Array{Float64}, wgx2::Array{Float64}, jac::Float64, x1gp::Int64, x2gp::Int64)

    res = 0.0 + 0.0im
    for l in 1:x2gp, k in 1:x1gp

        scale = wgx1[k] * wgx2[l] * jac

        for j in 1:9, i in 1:9
            res += @inbounds Ψ[i, k, l] * mat[i, j, k, l] * Φ[j, k, l] * scale
        end
    end

    return res
end



"""
    gauss_integrate(Ψ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, Φ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, mat::SubArray{ComplexF64, 4, Array{ComplexF64, 5}}, wgx1::Array{Float64}, wgx2::Array{Float64}, jac::Float64, x1gp::Int64, x2gp::Int64)
    
Computes the numerical integration required for finite elements via ∫f(x)dx ≈ ∑w_i f(x_i) in 3d.
"""
function gauss_integrate(Ψ::SubArray{ComplexF64, 4, Array{ComplexF64, 7}}, Φ::SubArray{ComplexF64, 4, Array{ComplexF64, 7}}, mat::Array{ComplexF64, 5}, wgx1::Array{Float64}, wgx2::Array{Float64}, wgx3::Array{Float64}, jac::Float64, x1gp::Int64, x2gp::Int64, x3gp::Int64)

    res = 0.0 + 0.0im
    for x3 in 1:x3gp, x2 in 1:x2gp, x1 in 1:x1gp

        scale = wgx1[x1] * wgx2[x2] * wgx3[x3] * jac

        for j in 1:9, i in 1:9
            res += @inbounds mat[i, j, x1, x2, x3] * scale * Ψ[i, x1, x2, x3] * Φ[j, x1, x2, x3]
        end
    end

    return res
end
