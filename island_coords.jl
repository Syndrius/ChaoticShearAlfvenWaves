
#understanding these bad bois based on mapping back, perhaps.

#we will use exactly Zhisongs form, except with r.
function tor_from_isl(κ, ab, φ, isl)

    if κ > 1
        #may change later!
        return 0.5, 0, 0, 0
    end

    α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ)/π * ab, κ))

    χ = -(2*isl.A*κ - isl.A)

    θ = α + φ/isl.q0

    #don't love the abs
    res = sqrt(abs(-(2*isl.q0^2)/ isl.qp *(χ - isl.A * cos(isl.m0 * α))))#+ isl.ψ0

    αm = mod(ab, 2π)

    #this is a good start, but this is arbitrary, not sure if this actually agrees with our alphabar definition.
    if αm < π/2
        ψ = -res + isl.ψ0
    elseif αm < 3π/2
        ψ = res + isl.ψ0
    else
        ψ = -res + isl.ψ0
    end

    r = sqrt(2*ψ)

    ζ = φ

    return r, θ, ζ, 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ)/π * ab, κ))

end



isl = ContIslandT(m0=3, n0=2, A=A, q0=3/2, qp=1.6, ψ0=0.125)



Nκ = 200
Nab = 80
Nφ = 30

κgrid = LinRange(0, 1.5, Nκ)
abgrid = LinRange(0, 2π, Nab+1)[1:end-1]
φgrid = LinRange(0, 2π, Nφ+1)[1:end-1]

rvals = zeros(Nκ, Nab, Nφ);
θvals = zeros(Nκ, Nab, Nφ);
ζvals = zeros(Nκ, Nab, Nφ);
αvals = zeros(Nκ, Nab, Nφ);

for (i, κ) in enumerate(κgrid), (j, ab) in enumerate(abgrid), (k, φ) in enumerate(φgrid)
    rvals[i, j, k], θvals[i, j, k], ζvals[i, j, k], αvals[i, j, k] = tor_from_isl(κ, ab, φ, isl)
end

contourf(abgrid, κgrid, rvals[:, :, 1] .- 0.5, levels=100, color=:turbo)
contourf(abgrid, κgrid, θvals[:, :, 1], levels=100, color=:turbo)
contourf(abgrid, κgrid, αvals[:, :, 10], levels=100, color=:turbo)
