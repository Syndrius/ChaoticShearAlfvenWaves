
"""
Function to compute ψ̄, used for plotting.

# Args
isl::ContIslandT - Struct storing island parameters.
χ::Array{Float64} - Energy values to consider.
sign::Int64 - Sign of the particles, ±1 for passing on either side, or 0 for trapped particles.
"""
function compute_ψ̄(isl::ContIslandT, χ::Array{Float64}, sign::Int64)
    κ = @. (isl.A - χ) / (2*isl.A)
    w = 4 * sqrt(isl.A  * isl.q0^2 /isl.qp)

    if sign == 0
        return @. 2*w/(isl.m0*π) * ((κ-1)*Elliptic.K(κ) + Elliptic.E(κ))
    else
        #we may one day want Ω, not sure though!
        #think it acts as the q-profile under new coords.
        #Ω = @. -sign * sqrt(isl.A * isl.qp / isl.q0^2) * π * sqrt(κ) / real(Elliptic.K(1/κ))
        return @. isl.ψ0 + sign * w / π * sqrt(κ) * real(Elliptic.E(1/κ))
    end
end

"""
Compute α in the original toroidal coordinates.

# Args
χ::Float64 - Energy value.
ᾱ::Float64 - α in new coordinates to be transformed.
isl::ContIslandT - Struct storing island parameters.
sign::Int64 - Sign of the particles, ±1 for passing on either side, or 0 for trapped particles.
"""
function compute_α(χ::Float64, ᾱ::Float64, isl::ContIslandT, sign::Int64)

    κ = (-χ + isl.A) / (2*isl.A)

    
    #could swap to if κ > 1 .
    if sign == 0
        #real seems bad.
        K = real.(Elliptic.K(κ))
        sn = real(sin(Elliptic.Jacobi.am(2*K /π * ᾱ, κ)))
        return 2/isl.m0 * asin.(sqrt.(κ)*sn)
    else
        K = real(Elliptic.K(1/κ))
        amval = Jacobi.am(isl.m0*K * ᾱ /π, 1/κ)
        return 2 / isl.m0 * amval
    end
end

"""
Compute ψ in the original toroidal coordinates.

# Args
χ::Float64 - Energy value.
α::Float64 - α in old coordinates.
ᾱ::Float64 - α in new coordinates.
isl::ContIslandT - Struct storing island parameters.
sign::Int64 - Sign of the particles, ±1 for passing on either side, or 0 for trapped particles.
"""
function compute_ψ(χ::Float64, α::Float64, ᾱ::Float64, isl::ContIslandT, sign::Int64)

    if sign == 0
        #would be good to know what the fuck this is doing.
        fac = 2*mod.( floor.((ᾱ .- π/2) / π), 2) .- 1
    else
        fac = sign
    end
    #order for cos - χ is different in trapped vs passing, abs takes care of this.
    #fac tells us which sign it should take.
    return fac * sqrt(2*isl.q0^2/isl.qp) * sqrt(abs(isl.A * cos(isl.m0*α) - χ)) + isl.ψ0 
end

        
"""
Compute ∇χ^2, the main term in the weakform.

# Args
ψ::Float64 - Radial coordinate.
α::Float64 - Helical angle.
met::MetT - metric.
isl::ContIslandT - Struct storing island parameters.
"""
function compute_∇χ2(ψ::Float64, α::Float64, met::MetT, isl::ContIslandT)

    return (isl.qp ^2 /isl.q0^4 * (ψ - isl.ψ0)^2 * met.gu[1, 1] 
            + isl.A^2 * isl.m0^2 *sin(α *isl.m0)^2 * (met.gu[2, 2] 
            + met.gu[3, 3] / isl.q0^2) 
            + 2*isl.qp / isl.q0^2 * (ψ - isl.ψ0) * isl.A * isl.m0 *sin(isl.m0*α) * met.gu[1, 2])
end


"""
Computes Ω, the local rotational transform.

# Args
χ::Float64 - Energy value.
isl::ContIslandT - Struct storing island parameters.
sign::Int64 - Sign of the particles, ±1 for passing on either side, or 0 for trapped particles.
"""
function compute_Ω(χ::Float64, isl::ContIslandT, sign::Int64)
    Ω0 = sqrt(isl.A * isl.qp / isl.q0^2)

    κ = (-χ + isl.A)/(2*isl.A)

    if sign == 0
        K = real.(Elliptic.K(κ))
        return -Ω0 * isl.m0*π / (2*K) 
    else
        return -sign * Ω0 * π * sqrt(κ)/Elliptic.K(1/κ)
    end
end


#inputs are cooked becuase we mostly want to use the normal island for this!
#probably possible to automate this with probs passed in, as then we could compute qp from dq
#and we can get q0 from m0/n0 or whatever.
function island_width(isl::IslandT, q0, qp)
    #with our form of island, the width is cooked, need some kind of solver for this!

    #for v small island, we can fix r(1-r) ≈ 0.25, and use normal width calculation,
    return 4 * sqrt(isl.A * 0.25 * q0^2 / qp)

end