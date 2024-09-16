
#handles the requirement to map between from toriodal coordinates
#into island coordinates

#needs a different name by golly.
function coords_isl_to_tor(κ, ᾱ, φ, isl::ContIslandT)
    #name of this is a bit confusing
    #but idea is that we pass in island coordinates
    #and we find equivalent toroidal coordinates.

    if κ > 1
        α = 2/isl.m0 * Elliptic.Jacobi.am(isl.m0 * Elliptic.K(1/κ) * ᾱ / π, 1/κ)
    else
        α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ))
    end

    χ = -1 * (2*isl.A * κ - isl.A)

    #lol wot even is r0...
    #taking abs is bold here, but only negatives are like e-20
    ψ = sqrt(abs(-2*isl.q0^2/isl.qp * (χ - isl.A*cos(isl.m0 * α)))) + isl.ψ0

    θ = α + φ/isl.q0
    #no fkn idea if this is actually flux or not lol.
    return ψ, θ, φ

end