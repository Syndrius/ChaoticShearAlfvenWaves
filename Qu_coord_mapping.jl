
#simplify this to only map a function directly from toroidal to island.

using Elliptic
using Plots; plotlyjs()
using MID
using Interpolations




isl = ContIslandT(m0=3, n0=2, qp=1.6, q0=3/2, A = 1e-3, ψ0 = 0.125)


function isl_to_toroidal(κ, ᾱ, φ, isl)
    #will this return r or flux??? who the fek knows.
    #problemo may be that κ is defined in terms of r^2.
    #so plotting it vs r makes no sense.
    #this is going to cause us issues big time.
    #κ = sqrt(κ)
    #κ = κ^2

    if κ > 1
        α = 2/isl.m0 * Elliptic.Jacobi.am(isl.m0 * Elliptic.K(1/κ) * ᾱ / π, 1/κ)
    else
        α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ))
    end

    χ = -1 * (2*isl.A * κ - isl.A)

    #lol wot even is r0...
    #taking abs is bold here, but only negatives are like e-20
    ψ = sqrt(abs(-2*isl.q0^2/isl.qp * (χ - isl.A*cos(isl.m0 * α)))) + isl.ψ0

    θ = α - φ/isl.q0

    return sqrt(ψ*2), θ, φ

end


function island_mode(r, θ, ζ, isl)

    r0 = sqrt(isl.ψ0*2)

    α = θ + ζ/isl.q0 #sign of this is questionable.
    χ = -isl.qp/(2*isl.q0^2) * (r^2/2 - r0^2/2)^2 + isl.A * cos(isl.m0*α)
    #χ = -isl.qp/(2*isl.q0^2) * (ψ - isl.r0)^2 + isl.A * cos(isl.m0*α)
    #κ = sqrt((-χ + isl.A) / (2*isl.A))
    κ = ((-χ + isl.A) / (2*isl.A))

    #outside the island!
    #don't seem to be having the numerical issues anymore,
    #I guess this is a just a coincidence since our grid is not hitting 0.5 anymore.
    if κ > 1
        #ᾱ = π/(2 * Elliptic.K(1/κ)) * Elliptic.F(isl.m0*α/2, 1/κ)
        #I guess this does have to be m0 not 2. not ideal...
        ᾱ = π/(isl.m0 * Elliptic.K(1/κ)) * Elliptic.F(isl.m0*α/2, 1/κ)
    else
        
        ᾱ = π/(2 * Elliptic.K(κ)) * Elliptic.F(asin(1/sqrt(κ)*sin(isl.m0*α/2)), κ)
    end

    #return κ, ᾱ, ζ

    m=2
    if κ > 1
        return 0
    end
    #return exp(-(κ-0.7)^2/0.01) * exp(1im * m * ᾱ)
    return (-4*(κ-0.5)^2+1)*exp(1im * m * ᾱ)
    #return (-4*(κ-0.5)^2+1)*cos(m * ᾱ)
end


Nr = 200
Nθ = 80
Nζ = 30

rgrid = LinRange(0, 1, Nr)
θgrid = LinRange(0, 2*π, Nθ+1)[1:end-1]
ζgrid = LinRange(0, 2π, Nζ+1)[1:end-1]

ϕ_tor = zeros(ComplexF64, Nr, Nθ, Nζ);

for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)
    ϕ_tor[i, j, k] = island_mode(r, θ, ζ, isl)
end

contourf(θgrid, rgrid, real.(ϕ_tor[:, :, 1]), levels=100, color=:turbo)


tor_itp = interpolate((rgrid, θgrid, ζgrid), ϕ_tor, Gridded(Linear(Periodic())));

tor_ext = extrapolate(tor_itp, Periodic());

Nκ2 = 200
Nᾱ2 = 80
Nφ2 = 30
#κgrid2 = LinRange(0, 8, Nκ2);
κgrid2 = LinRange(0, 1.5, Nκ2);
ᾱgrid2 = LinRange(0, 2π, Nᾱ2+1)[1:end-1];
φgrid2 = LinRange(0, 2π, Nφ2+1)[1:end-1];

ϕ_isl2 = zeros(ComplexF64, Nκ2, Nᾱ2, Nφ2);


for (i, κ) in enumerate(κgrid2), (j, ᾱ) in enumerate(ᾱgrid2), (k, φ) in enumerate(φgrid2)

    r, θ, ζ = isl_to_toroidal(κ, ᾱ, φ, isl)
    
    ϕ_isl2[i, j, k] = tor_ext(r, θ, ζ)

end

#does work, but we may be better of with Axel's coordinates.
contourf(ᾱgrid2, κgrid2, real.(ϕ_isl2[:, :, 1]), levels=100, color=:turbo)