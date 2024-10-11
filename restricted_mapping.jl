
#testing some of our mapping, but with restricted domains.



using Elliptic
using Plots; plotlyjs()
using MID
using Interpolations
using FFTW
using Printf




isl = IslandT(m0=2, n0=1, qp=1.6, q0=3/2, A = 1e-3, r0 = 0.125)


function isl_to_tor(κ, ᾱ, φ, isl)
    #will this return r or flux??? who the fek knows.
    #problemo may be that κ is defined in terms of r^2.
    #so plotting it vs r makes no sense.
    #this is going to cause us issues big time.
    #κ = sqrt(κ)
    #κ = κ^2
    αmod = mod(ᾱ, 2π) #probably not needed as we are passing in (0, 2π).

    if κ > 1
        α = 2/isl.m0 * Elliptic.Jacobi.am(isl.m0 * Elliptic.K(1/κ) * ᾱ / π, 1/κ)
    else

        α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ))
    end

    χ = -1 * (2*isl.A * κ - isl.A)

    #lol wot even is r0...
    #taking abs is bold here, but only negatives are like e-20

    #maybe we need to mod α, unsure if we need to mod alpha or m0 alpha though!
    res = sqrt(abs(-2*isl.q0^2/isl.qp * (χ - isl.A*cos(isl.m0 * α))))

    
    #r0 = sqrt(2*isl.ψ0)
    #not sure if ψ0 and ψ should be converted together or individiually.
    if αmod < π/2
        ψ = +res + isl.ψ0
        #r = sqrt(2*res) + r0
    elseif αmod < 3π/2
        ψ = -res + isl.ψ0
        #r = -sqrt(2*res) + r0
    else
        ψ = +res + isl.ψ0
        #r = sqrt(2*res) + r0
    end

    θ = α + φ/isl.q0
    #display(α)
    try 
        r = sqrt(ψ*2)
    catch 
        #seems to only occur when the island is to big??? v odd.
        display("Sqrt no good")
        display(α)
    end
    return r, mod(θ, 2π), mod(φ, 2π)
    #return r, mod(θ, 2π), mod(φ, 2π)

end

function tor_to_isl(r, θ, ζ, isl)
    r0 = sqrt(isl.r0*2)

    #taking mod here cooks this apparently...
    #maybe a hint of what is going on??? Perhaps the elliptic stuff needs (-pi, pi)?
    #unsure why that would be????
    #α = mod(θ - ζ/isl.q0, 2π) #sign of this is questionable.
    α = θ - ζ/isl.q0 #sign of this is questionable.
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
        #ᾱ=0
    else

        #if α > π
        #    α = α - π
        #end
        amod = mod(α * isl.m0, 2π)
        #This set up gives us what we want, but obviously the abs is less than ideal
        if r > r0
            if mod(α * isl.m0, 2π) < π
        
                #ᾱ = π/(2 * Elliptic.K(κ)) * Elliptic.F(-asin(1/sqrt(κ)*sin(isl.m0*α/2))+π, κ)
                ᾱ = π/(2 * Elliptic.K(κ)) * Elliptic.F(asin(1/sqrt(κ)*abs(sin(isl.m0*α/2))), κ)# + π/2
                #ᾱ = π/(2 * Elliptic.K(κ)) * Elliptic.F(-asin(1/sqrt(κ)*sin(isl.m0*amod/2)), κ)# + π/2
                #ᾱ = π/(2 * Elliptic.K(κ)) * Elliptic.F(asin(1/sqrt(κ)*sin(isl.m0*α/2)), κ)
            else
                ᾱ = π/(2 * Elliptic.K(κ)) * Elliptic.F(-asin(1/sqrt(κ)*abs(sin(isl.m0*α/2)))+2π, κ)
                #ᾱ=0
            end

        else
            if mod(α* isl.m0, 2π) < π
                ᾱ = π/(2 * Elliptic.K(κ)) * Elliptic.F(-asin(1/sqrt(κ)*abs(sin(isl.m0*α/2))) + π, κ)# - π/2
                #ᾱ=0
            else
                ᾱ = π/(2 * Elliptic.K(κ)) * Elliptic.F(asin(1/sqrt(κ)*abs(sin(isl.m0*α/2)))+π, κ)# - π/2
                #ᾱ=0
            end
        end
    end

    #return κ, ᾱ, ζ

    ᾱ = mod(ᾱ, 2π)

    return κ, ᾱ, ζ
end


function inside_island_mode(κ, ᾱ, φ)

    m=2
    n = 1
    if κ > 1
        return 0
    end

    if 2π/2 < mod(ᾱ, 2π) #< 2π #useful for resricting the mapping.
        return 0
    end

    #return exp(-(κ-0.7)^2/0.01) * exp(1im * m * ᾱ)
    return (-4*(κ-0.5)^2+1)*exp(1im * m * ᾱ + 1im * n * φ/3)
    #return (-4*(κ-0.5)^2+1)*cos(m * ᾱ)
end

function island_mode(r, θ, ζ, isl)

    #if π/2 < θ < 3π/2
    #    return 0
    #end

    if ζ > 10π/3
        return 0
    end

    κ, ᾱ, φ = tor_to_isl(r, θ, ζ, isl)

    m=2
    n = 3
    if κ > 1
        return 0
    end

    if 5π/2 < mod(ᾱ, 2π)# < π #useful for resricting the mapping.
        return 0
    end
    #return ᾱ
    #return exp(-(κ-0.5)^2/0.01) * exp(1im * m * ᾱ  + 1im * n * φ)
    return (-4*(κ-0.5)^2+1)*exp(1im * m * ᾱ - 1im * n * φ/3)
    #return (-4*(κ-0.5)^2+1)*cos(m * ᾱ)
end

"""
Nκ = 200
Nᾱ = 80
Nφ = 30
#κgrid = LinRange(0, 8, Nκ);
κgrid = LinRange(0, 1.5, Nκ);
ᾱgrid = LinRange(-π, π, Nᾱ)#+1)[1:end-1];
φgrid = LinRange(-π, π, Nφ)#+1)[1:end-1];


ϕ_isl = zeros(ComplexF64, length(κgrid), length(ᾱgrid), length(φgrid));




for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, φ) in enumerate(φgrid)
    ϕ_isl[i, j, k] = inside_island_mode(κ, ᾱ, φ)
end

contourf(ᾱgrid, κgrid, real.(ϕ_isl[:, :, 10]), color=:turbo, levels=100)
""";


Nr = 200
Nθ = 80
Nζ = 30

rgrid = LinRange(0, 1, Nr)
θgrid = LinRange(0, 2π, Nθ)#+1)[1:end-1]
ζgrid = LinRange(0, 2π, Nζ)#+1)[1:end-1]

ϕ_tor = zeros(ComplexF64, Nr, Nθ, Nζ);

κvals = zeros(Nr, Nθ, Nζ);
ᾱvals = zeros(Nr, Nθ, Nζ);
αvals_tor = zeros(Nr, Nθ, Nζ);

for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)
    αvals_tor[i, j, k] = θ - ζ/isl.q0
    κvals[i, j, k], ᾱvals[i, j, k], _ = tor_to_isl(r, θ, ζ, isl)
    ϕ_tor[i, j, k] = island_mode(r, θ, ζ, isl)
end


contourf(θgrid, rgrid, ᾱvals[:, :, 10], levels=100, color=:turbo)
contourf(θgrid, rgrid, mod.(2 .* ᾱvals[:, :, 1], 2π), levels=100, color=:turbo)

contourf(θgrid, rgrid, mod.(5 .* αvals_tor[:, :, 1], 2π), levels=100, color=:turbo)
contourf(θgrid, rgrid,  αvals_tor[:, :, 28], levels=100, color=:turbo)

#ϕ_tor = periodify(ϕ_tor_a, Nr, Nθ, Nζ);


contourf(θgrid, rgrid, real.(ϕ_tor[:, :, 5]), levels=100, color=:turbo)


itp_tor = interpolate((rgrid, θgrid, collect(ζgrid)), ϕ_tor, Gridded(Linear(Periodic())));

ext_tor = extrapolate(itp_tor, Periodic());

#now lets map this bad boi back!
Nκ2 = 200
Nᾱ2 = 100
Nφ2 = 30
#κgrid2 = LinRange(0, 8, Nκ2);
κgrid2 = LinRange(0, 1.5, Nκ2);
ᾱgrid2 = LinRange(0, 2π, Nᾱ2)#+1)[1:end-1];
φgrid2 = LinRange(0, 2π, Nφ2)#+1)[1:end-1];

ϕ_isl2 = zeros(ComplexF64, Nκ2, Nᾱ2, Nφ2);

θvals = zeros(Nκ2, Nᾱ2, Nφ2);
rvals = zeros(Nκ2, Nᾱ2, Nφ2);
ζvals = zeros(Nκ2, Nᾱ2, Nφ2);
αvals = zeros(Nκ2, Nᾱ2, Nφ2);


for (i, κ) in enumerate(κgrid2), (j, ᾱ) in enumerate(ᾱgrid2), (k, φ) in enumerate(φgrid2)

    r, θ, ζ = isl_to_tor(κ, ᾱ, φ, isl)

    rvals[i, j, k] = r
    θvals[i, j, k] = θ
    ζvals[i, j, k] = ζ
    αvals[i, j, k] = θ - ζ/isl.q0

    ϕ_isl2[i, j, k] = ext_tor(r, θ, ζ)

end

#so our island modes defintatly have a problem, they are not periodic in theta, which they definatly should be.

#(2, 1) mode may signify some other problemos, hard to be sure though!
#hard to tell what is going on tbh.
#including the periodicity does seem to be a big improvement.

contourf(ᾱgrid2, κgrid2, real.(ϕ_isl2[:, :, 1]), levels=100, color=:turbo)
contourf(ᾱgrid2, κgrid2, real.(θvals[:, :, 7]), levels=100, color=:turbo)
contourf(ᾱgrid2, κgrid2, real.(αvals[:, :, 7]), levels=100, color=:turbo)
contourf(ᾱgrid2, κgrid2, real.(rvals[:, :, 1] .- 0.5), levels=100, color=:turbo)


contourf(θgrid, rgrid, κvals[:, :, 1], levels=100, color=:turbo)
contourf(θgrid, rgrid, ᾱvals[:, :, 7], levels=100, color=:turbo)
contourf(θgrid, rgrid, mod.(αvals_tor[:, :, 8], 2π), levels=100, color=:turbo)
contourf(θgrid, rgrid, real.(ϕ_tor[:, :, 3]), levels=100, color=:turbo)

abgrid = LinRange(0, 2π, 101)[1:end-1]

κval = 0.7

plot(abgrid, @. 1/sqrt(κval)*sin(isl.m0*abgrid))


ϕ_islfft = fft(ϕ_isl2, [2, 3]);

p = plot()
for i in 1:Nᾱ, j in 1:Nφ

    if i > Nᾱ/2
        m = i - Nᾱ - 1
    else
        m = i-1
    end
    if j > Nφ/2
        n = j - Nφ - 1
    else
        n = j-1
    end

    plot!(κgrid, real.(ϕ_islfft[:, i, j]), label=@sprintf("(%s, %s)", m, n))
    #plot!(κgrid, imag.(ϕ_islfft[:, i, j]), label=@sprintf("(%s, %s)", m, n))
    #plot!(κgrid, real.(ϕisl[:, i, 1]))
end
#this may just be total garbage who knows.
display(p)





function res_ab(κ, θ, ζ, isl)

    α = θ - ζ/isl.q0

    α = θ - ζ/isl.q0 #sign of this is questionable.

    #outside the island!
    #don't seem to be having the numerical issues anymore,
    #I guess this is a just a coincidence since our grid is not hitting 0.5 anymore.
    if κ > 1
        #ᾱ = π/(2 * Elliptic.K(1/κ)) * Elliptic.F(isl.m0*α/2, 1/κ)
        #I guess this does have to be m0 not 2. not ideal...
        return 0
    else

        if abs((1/sqrt(κ)*sin(isl.m0*α/2))) > 1
            return 0
        end

        #if α > π
        #    α = α - π
        #end

        if mod(isl.m0*α, π) > π/2
        
            ᾱ = π/(2 * Elliptic.K(κ)) * Elliptic.F(-asin(1/sqrt(κ)*sin(isl.m0*α/2)), κ)
            #ᾱ = π/(2 * Elliptic.K(κ)) * Elliptic.F(asin(1/sqrt(κ)*sin(isl.m0*α/2)), κ)
        elseif mod(isl.m0*α, π) < 3π/2
            #ᾱ = π/(2 * Elliptic.K(κ)) * Elliptic.F(-asin(1/sqrt(κ)*sin(isl.m0*α/2))+π, κ)
            ᾱ = π/(2 * Elliptic.K(κ)) * Elliptic.F(asin(1/sqrt(κ)*sin(isl.m0*α/2)), κ) #+ π
        else
            ᾱ = π/(2 * Elliptic.K(κ)) * Elliptic.F(-asin(1/sqrt(κ)*sin(isl.m0*α/2)), κ)
        end
    end

    return ᾱ
end

function res_a(κ, ᾱ, φ, isl)
    if κ > 1
        α = 2/isl.m0 * Elliptic.Jacobi.am(isl.m0 * Elliptic.K(1/κ) * ᾱ / π, 1/κ)
    else
        α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ))
    end

    χ = -1 * (2*isl.A * κ - isl.A)



    #lol wot even is r0...
    #taking abs is bold here, but only negatives are like e-20
    res = sqrt(abs(-2*isl.q0^2/isl.qp * (χ - isl.A*cos(isl.m0 * α))))

    αmod = mod(ᾱ, 2π) #probably not needed as we are passing in (0, 2π).
    #r0 = sqrt(2*isl.ψ0)
    #not sure if ψ0 and ψ should be converted together or individiually.
    if αmod < π/2
        ψ = -res + isl.ψ0
        #r = sqrt(2*res) + r0
    elseif αmod < 3π/2
        ψ = res + isl.ψ0
        #r = -sqrt(2*res) + r0
    else
        ψ = -res + isl.ψ0
        #r = sqrt(2*res) + r0
    end

    θ = α + φ/isl.q0
    #no fkn idea if this is actually flux or not lol.
    #display(res)
    #display((κ,  αmod, α))
    #if κ == 0.0
    #    r = 0.5
    #else
        #r = sqrt(2π)
    #end
    return α
end


κval = 0.9999
ζval = 0.0
θlist = LinRange(0, 2π, 100)
ablist = LinRange(0, 2π, 100)
alist = zeros(100)
ablist = zeros(100)


for (i, abval) in enumerate(ablist)


    alist[i] = res_a(κval, abval, ζval, isl)
end

for (i, θval) in enumerate(θlist)


    ablist[i] = res_ab(κval, θval, ζval, isl)
end

#plot(θlist, alist)
plot(θlist, ablist)
plot(θlist, mod.(isl.m0*θlist, π))

#mayhap it would be better to work with -π to π, just for arcsin's sake!
α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ))

plot(ablist, @. mod(2 * asin(sqrt(κval)*Elliptic.Jacobi.sn(2*Elliptic.K(κval) * ablist / π, κval)), 2π))


#need to do this mapping wiht m=1 and m=-1 to see what changes. Would be good to understand why we get both/ what is happening when they are the same/different.

2 * sqrt(1e-4 * 1.6)