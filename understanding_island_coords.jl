
using MID
using Plots; plotlyjs()
using Elliptic


A = 1e-4

isl = ContIslandT(m0=3, n0=2, A=A, q0=3/2, qp=1.6, ψ0=0.125)#srs wot is psi 0


geo = GeoParamsT(R0=1000.0)

#start with this???
pmd = init_sm_grid(start=-12, count=26, incr=1)
tmd = init_sm_grid(start=-4, count=5, incr=2)
#tmd = MID.ModeDataT(start=-8, count=10, incr=2)

#χlist = LinRange(-A+A*0.05, A-A*0.05, 50) #this is fked, think we need to cluster near the spratrix.

#just try what Zhisong did?
χlist1 =  1 .- 0.01 * exp.( -1 .* collect(LinRange(0, 12, 11)))
χlist2 = LinRange(0, sqrt(χlist1[1]), 191)[1:end-1] .^2 .+ 1e-5
#ok this matches Zhisong, again no fkn idea what the hell this list is.
χlist = A .- vcat(χlist2, χlist1) .* 2 .* A

ω2list = island_continuum(χlist, pmd, tmd, geo, isl, 0);

#width = 4 * sqrt(isl.A * isl.q0^2/isl.qp)
#ψ_isl = 2 * width / (π * isl.m0)
#ψ̄m = ψ_isl
ψ̄list = MID.Continuum.compute_ψ̄(isl, χlist, 0)

r = repeat(sqrt.(2 .* ψ̄list), 1,  pmd.count * tmd.count)

#we should normalise the x-axis, it is going from centre of island to edge I think.
scatter(r, sqrt.(abs.(ω2list .* geo.R0^2)), legend=false, markersize=0.1, ylimits=(0.0, 0.1))#.2, 0.6))



#lets consider Axel's in between coords

r0 = 0.5
w = 0.2 #defined somewhat arbitrarily.
m0 = 3
n0 = 2

function bent_axel_kappa(r, θ, ζ)

    α = 1/2 * (m0*θ + n0*ζ)
    return sqrt(1/w^2 * (r^2-r0^2)^2 + sin(α)^2)
end

#absolutely no understanding of this.
function bent_axel_beta(r, θ, ζ)

    α = 1/2 * (m0*θ + n0*ζ)

    #so this makes a huge difference!
    return atan(w*sin(α), (r^2-r0^2))
    #return atan(w*sin(α) / (r^2-r0^2))
end

function straight_axel_beta(κ, β)
    if κ < 1
        return π / (2 * Elliptic.K(κ)) * Elliptic.F(β, κ)
    else
        return 0
    end
end

function alpha_bar(r, θ, ζ)
    #this should be identical to βstraight
    
    A = 0.003 #subject to change obvs.
    q0 = 3/2
    α = θ - ζ / q0
    qp = 4
    χ = -qp / (2 * q0^2) * (r^2/2 - r0^2/2)^2 + A*cos(m0 * α)
    #think Axels kappa is like a prenormalised χ.
    κ = (-χ + A) / (2 * A)

    if κ < 1

        #return 1/sqrt(κ) * sin(m0 * α / 2)
        #return mod(m0 * α/2, 2π)

        ab = π / (2 * Elliptic.K(κ)) * Elliptic.F(asin(1/sqrt(κ) * sin(m0 * α / 2)), κ)

        #return ab
        #return mod(ab, 2π)
        if r > r0
            ab = π / (2 * Elliptic.K(κ)) * Elliptic.F(-asin(1/sqrt(κ) * sin(m0 * α / 2)) + π, κ)
            #return ab #- 3π/2
            return mod(ab + 0*π, 2π)
        else
            #think this is it, note that we are just changin asin to -asin + π.
            #we can swap them around, just dictates if ᾱ=0 is inner or outer part of torus.
            #not sure if one way would be better than the other??
            ab = π / (2 * Elliptic.K(κ)) * Elliptic.F(asin(1/sqrt(κ) * sin(m0 * α / 2)), κ)
            #return ab #+ π
            return mod(ab+0*3π/2, 2π)
        end

        #return π / (2 * Elliptic.K(κ)) * Elliptic.F(asin(1/sqrt(κ) * sin(m0 * α / 2)), κ)
        return π / (2 * Elliptic.K(κ)) * Elliptic.F(asin(1/sqrt(κ) * sin(m0 * α / 2)), κ)
        #return π / (2 * Elliptic.K(κ)) * 1/sqrt(κ) * sin(m0 * α / 2)
    else 
        return 0
    end
end



Nr = 100
Nθ = 80
Nζ = 10

rgrid = LinRange(0, 1, Nr)
θgrid = LinRange(0, 2π, Nθ+1)[1:end-1]
ζgrid = LinRange(0, 2π, Nζ+1)[1:end-1]

κbent = zeros(Nr, Nθ, Nζ);
βbent = zeros(Nr, Nθ, Nζ);
βstraight = zeros(Nr, Nθ, Nζ);
ᾱ = zeros(Nr, Nθ, Nζ);
α = zeros(Nr, Nθ, Nζ);
for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)
    κ = bent_axel_kappa(r, θ, ζ)
    β = bent_axel_beta(r, θ, ζ)
    κbent[i, j, k] = κ
    βbent[i, j, k] = β
    βstraight[i, j, k] = straight_axel_beta(κ, β)
    ᾱ[i, j, k] = alpha_bar(r, θ, ζ)
    α[i, j, k] = mod(3*θ - 2*ζ, 2π)
    #α[i, j, k] = mod(θ + 2*ζ/3, 2π)

end

contourf(θgrid, rgrid, ᾱ[:, :, 1], levels=100)

#island is non-symmetrical.... interesti

#we may need to start mapping these to a torus to undertsand, there is probably some significance for inner vs outer...

contourf(θgrid, rgrid, κbent[:, :, 1], levels=100)
contourf(θgrid, rgrid, βbent[:, :, 1], levels=200)

#this makes much more sense as far as an island is concerned
contourf(θgrid, rgrid, βstraight[:, :, 1], levels=100)

contourf(θgrid, rgrid, α[:, :, 1], levels=100)


scatter(Zgrid, Rgrid, tor_func, levels=200)

function convert_to_torus(rgrid, θgrid, ζgrid, func)
    #this could be done with interpolation, but seems like it shouldn't be needed???
    #defo is big rip.
    #don't have inverse transformation??? going to be super cooked!

    #surely interpolation with root finding is not the best way to do this???
    #func is assumed 3d

    #going to try a real sloppy method...
    #Rgrid = zeros(length(rgrid), length(θgrid))
    #Zgrid = zeros(length(rgrid), length(θgrid))

    Rgrid = []
    Zgrid = []

    tor_func = []

    #may ignore zeta for now.
    for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid)

        Rval, _, Zval = radial_to_plot(r, θ, 1, 1000)
        push!(Rgrid, Rval)
        push!(Zgrid, Zval)
        push!(tor_func, func[i, j, 1])


    end

    return Rgrid, Zgrid, tor_func


end


Rgrid, Zgrid, tor_func = convert_to_torus(rgrid, θgrid, ζgrid, βbent)

function radial_to_plot(r, θ, ζ, R0)
    #now we are actually using flux.
    #display(ψ)

    Δ = @. r^2/(R0 * 8) - 1/(R0*8)
    Δp = @. r / (R0 * 4)
    η = @. (r/R0 +Δp)/2

    R = @. R0 + r*cos(θ) - Δ + r * η * (cos(2*θ) - 1)
    Z = @. r * sin(θ) + r * η * sin(2*θ)
    return [R, -ζ, Z]
end