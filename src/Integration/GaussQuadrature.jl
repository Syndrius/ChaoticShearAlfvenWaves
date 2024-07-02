

"""
    gauss_integrate(Ψ::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, Φ::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, mat::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, wg::Array{Float64}, jac::Float64, ngp::Int64)

Computes the numerical integration required for finite elements via ∫f(x)dx ≈ ∑w_i f(x_i) in 1d.
"""
function gauss_integrate(Ψ::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, Φ::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, mat::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, wg::Array{Float64}, jac::Float64, ngp::Int64)

    #this function is taking up a lot of time when nm, nn get large, but I think that is more to do with the size of the loops that call this function tbh.

    res = 0.0 + 0.0im
    for k in 1:ngp
    
        scale = wg[k] * jac
        for j in 1:9, i in 1:9
            #significantly faster to have * jac here not later!
            res += @inbounds Ψ[i, k] * mat[i, j, k] * Φ[j, k] * scale

        end
    end
    #res *= jac
    return res #* dr / 2 
end



"""
    gauss_integrate(Ψ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, Φ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, mat::SubArray{ComplexF64, 4, Array{ComplexF64, 5}}, wgr::Array{Float64}, wgθ::Array{Float64}, jac::Float64, rgp::Int64, θgp::Int64)
    
Computes the numerical integration required for finite elements via ∫f(x)dx ≈ ∑w_i f(x_i) in 2d.
"""
function gauss_integrate(Ψ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, Φ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, mat::SubArray{ComplexF64, 4, Array{ComplexF64, 5}}, wgr::Array{Float64}, wgθ::Array{Float64}, jac::Float64, rgp::Int64, θgp::Int64)

    res = 0.0 + 0.0im
    for l in 1:θgp, k in 1:rgp

        scale = wgr[k] * wgθ[l] * jac

        for j in 1:9, i in 1:9
            res += @inbounds Ψ[i, k, l] * mat[i, j, k, l] * Φ[j, k, l] * scale
        end
    end

    return res
end



"""
    gauss_integrate(Ψ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, Φ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, mat::SubArray{ComplexF64, 4, Array{ComplexF64, 5}}, wgr::Array{Float64}, wgθ::Array{Float64}, jac::Float64, rgp::Int64, θgp::Int64)
    
Computes the numerical integration required for finite elements via ∫f(x)dx ≈ ∑w_i f(x_i) in 2d.
"""
function gauss_integrate(Ψ::SubArray{ComplexF64, 4, Array{ComplexF64, 7}}, Φ::SubArray{ComplexF64, 4, Array{ComplexF64, 7}}, mat::Array{ComplexF64, 5}, wgr::Array{Float64}, wgθ::Array{Float64}, wgζ::Array{Float64}, jac::Float64, rgp::Int64, θgp::Int64, ζgp::Int64)

    res = 0.0 + 0.0im
    for ζ in 1:ζgp, θ in 1:θgp, r in 1:rgp

        scale = wgr[r] * wgθ[θ] * wgζ[ζ] * jac

        for j in 1:9, i in 1:9
            res += @inbounds Ψ[i, r, θ, ζ] * mat[i, j, r, θ, ζ] * Φ[j, r, θ, ζ] * scale
        end
    end

    return res
end
