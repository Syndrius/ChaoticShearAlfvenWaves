
#Currently this file serves no purpose!



"""
Computes the numerical integration required for finite elements via ∫f(x)dx ≈ ∑w_i f(x_i).
#this function contracts the matrix with the basis we have made and integrates using gauss quadrature.
#always 9x9 now.
#these are extra stupid types as we are taking views, we have just copied the error message, this may be a bad idea!

# Args These are stupid af and probably subject to change.
"""
function gauss_integrate(Ψ::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, Φ::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, mat::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, wg::Array{Float64}, jac::Float64, ngp::Int64)#::ComplexF64

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


#need to know the types of these bad bois!
#this version if for 2d fem.
function gauss_integrate(Ψ, Φ, mat, wgr, wgθ, jac, rgp, θgp)

    res = 0.0 + 0.0im
    for l in 1:θgp, k in 1:rgp

        scale = wgr[k] * wgθ[l] * jac

        for j in 1:9, i in 1:9
            res += @inbounds Ψ[i, k, l] * mat[i, j, k, l] * Φ[j, k, l] * scale
        end
    end

    return res
end


#think this is useless now!
function gauss_integrate_for_big(test_vec::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, trial_vec::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, mat::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, wg::Array{Float64}, jac::Float64, ngp::Int64)#::ComplexF64

    #this function is taking up a lot of time when nm, nn get large, but I think that is more to do with the size of the loops that call this function tbh.

    res = 0.0 + 0.0im
    for k in 1:ngp
    
        scale = wg[k] * jac
        for j in 1:10

            for i in 1:10
                #significantly faster to have * jac here not later!
                res += @inbounds test_vec[i, k] * mat[i, j, k] * trial_vec[j, k] * scale
            #display(res)
            end
        end
    end
    #res *= jac
    return res #* dr / 2 
end

#=
#not used anymore, checked if it was faster, but it is not, kinda surprised tbh!
function combined_integrate(resI, resW, test_vec, trial_vec, I, W, wg, jac, ngp, left_dim, right_dim)#::Tuple{ComplexF64, ComplexF64}
    #this one is a wee bit slower for some reason, has more gc time..
    resI = 0.0 + 0.0im
    resW = 0.0 + 0.0im
    for i in 1:left_dim

        for j in 1:right_dim

            for k in 1:ngp

                resI += @inbounds @views test_vec[i, k] * I[i, j, k] * trial_vec[j, k] * wg[k] * jac
                resW += @inbounds @views test_vec[i, k] * W[i, j, k] * trial_vec[j, k] * wg[k] * jac
            #display(res)
            end
        end
    end
    #return (resI * dr / 2, resW * dr / 2)

end
=#

