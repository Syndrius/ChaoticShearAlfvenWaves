

function κ(r, θ, ζ, isl)

    r0 = sqrt(isl.ψ0*2)

    #taking mod here cooks this apparently...
    #maybe a hint of what is going on??? Perhaps the elliptic stuff needs (-pi, pi)?
    #unsure why that would be????
    #α = mod(θ - ζ/isl.q0, 2π) #sign of this is questionable.
    α = θ - ζ/isl.q0 #sign of this is questionable.
    #χ = -isl.qp/(2*isl.q0^2) * (r^2/2 - r0^2/2)^2 + isl.A * cos(isl.m0*α)
    χ = -isl.qp/(2*isl.q0^2) * (r - r0)^2 + isl.A * cos(isl.m0*α)
    #χ = -isl.qp/(2*isl.q0^2) * (ψ - isl.r0)^2 + isl.A * cos(isl.m0*α)
    #κ = sqrt((-χ + isl.A) / (2*isl.A))
    κ = ((-χ + isl.A) / (2*isl.A))

    return κ
end


function κ_tor(κlist)

    #will ignore zeta for now
    rvals = zeros(length(κlist))
    θvals = zeros(length(κlist))
    #κvals = zeros(length(κlist), length(κlist))
    for κ in κlist
        α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ))

    end

    return rvals, θvals#, ζvals


end

        


