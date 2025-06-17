
"""
    function compute_I!(I::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, B::BFieldT, met::MetT, n::Float64, flr::FLRT, D::Array{Float64, 2}, F::Array{Float64})

Computes the I matrix for the weak form at a single coordinate. 
"""
function compute_I!(I::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, B::BFieldT, met::MetT, n::Float64, flr::FLRT, D::Array{Float64, 2}, F::Array{Float64})


    compute_F!(B, met, F)

    @inbounds for j=1:9, i=1:9
        #computes the non-ideal effects
        #assumes that Te ≈ Ti
        #then ρ_s = ρ_i
        #δ is artifical resitivity
        #ρ_i is ion gyro radius
        #δ_e is electron resistivity
        I[i, j] = -F[i] * F[j] * met.J[1] * n / B.mag_B[1]^2*(1.0im * flr.δ + 3/4 * flr.ρ_i^2 + (1-flr.δ_e*1im) * flr.ρ_i^2)
        I[i, j] = 0.0
    end
    
    #this is the ideal term, it is computed second to ensure I is zero everywhere first.
    I[1:3, 1:3] += D .* (n * met.J[1] / B.mag_B[1]^2)


end


"""
    compute_F!(B::BFieldT, met::MetT, F::Array{Float64})

Computes the F array, used for the non-ideal part of I.
"""
function compute_F!(B::BFieldT, met::MetT, F::Array{Float64})

    #TODO
    #need description of this, when we return to damping.

    F .= 0.0

    @inbounds for j=1:3, i=1:3
        F[j] += 1/met.J[1] * ( met.dJ[i]*(met.gu[i, j] - B.b[i]*B.b[j])
            + met.J[1]*(met.dgu[i, j, i] - B.b[i]*B.db[j, i] - B.db[i, i] * B.b[j]))
    end

    #ϕ_rr, i=r, j=r
    F[4] = met.gu[1, 1]-B.b[1]*B.b[1]
    #ϕ_rθ, i=r, j=θ + i=θ, j=r
    F[5] = met.gu[1, 2]-B.b[1]*B.b[2] + met.gu[2, 1]-B.b[2]*B.b[1]
    #ϕ_rζ, i=r, j=ζ + i=ζ, j=r
    F[6] = met.gu[1, 3]-B.b[1]*B.b[3] + met.gu[3, 1]-B.b[3]*B.b[1]
    #ϕ_θθ, i=θ, j=θ
    F[7] = met.gu[2, 2]-B.b[2]*B.b[2] 
    #ϕ_θζ, i=θ, j=ζ + i=ζ, j=θ
    F[8] = met.gu[2, 3]-B.b[2]*B.b[3] + met.gu[3, 2]-B.b[3]*B.b[2]
    #ϕ_ζζ, i=ζ, j=ζ
    F[9] = met.gu[3, 3]-B.b[3]*B.b[3] 

end
