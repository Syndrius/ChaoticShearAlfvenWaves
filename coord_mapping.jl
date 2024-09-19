
#lets try and do this properly
#i.e. fk off with psi. Hopefully that is the problemo..
using Elliptic
using Plots; plotlyjs()
using MID
using Interpolations

#this may have to be considered at at adequate standard...

#island type used for the island coordinate q-profile.
@kwdef struct QIslandT
    m0 :: Int64 
    n0 :: Int64
    A :: Float64
    q0 :: Float64
    qp :: Float64
    r0 :: Float64 #ideally change this to r if possible not sure if we can...
end

function inside_island_mode(κ, ᾱ, φ)

    m=2
    if κ > 1
        return 0
    end
    #return exp(-(κ-0.7)^2/0.01) * exp(1im * m * ᾱ)
    return (-4*(κ-0.5)^2+1)*exp(1im * m * ᾱ)
    #return (-4*(κ-0.5)^2+1)*cos(m * ᾱ)
end


function outside_island_mode(κ, ᾱ, φ)
    m=2
    if κ < 1
        return 0
    end
    #assumes κ goes to 8.
    return exp(-(κ-6.7)^2/0.01) * exp(1im * m * ᾱ)
    #return exp(-(κ-6.7)^2/0.01) * cos(m * ᾱ)
end

function create_island_mode(κgrid, ᾱgrid, φgrid, inside)

    

    ϕ_isl = zeros(ComplexF64, length(κgrid), length(ᾱgrid), length(φgrid))

    if inside
        


        for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, φ) in enumerate(φgrid)
            ϕ_isl[i, j, k] = inside_island_mode(κ, ᾱ, φ)
        end
    else
        #κgrid = LinRange(0, 5, Nκ)


        for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, φ) in enumerate(φgrid)
            ϕ_isl[i, j, k] = outside_island_mode(κ, ᾱ, φ)
        end

    end
    return ϕ_isl

end


#finds the equivalent (κ, ᾱ, φ) for a given (r, θ, ζ)
function toroidal_to_isl(ψ, θ, ζ, isl)

    #so just pretending r is ψ removed the symmetry problemos.
    #may just have to work with ψ...
    α = θ - ζ/isl.q0
    χ = -isl.qp/(2*isl.q0^2) * (ψ - isl.r0^2/2)^2 + isl.A * cos(isl.m0*α)
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

    return κ, ᾱ, ζ

end

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
    ψ = sqrt(abs(-2*isl.q0^2/isl.qp * (χ - isl.A*cos(isl.m0 * α)))) + isl.r0

    θ = α + φ/isl.q0

    return ψ, θ, φ

end


#first we create our island mode
isl = QIslandT(m0=3, n0=2, A=1e-3, q0=3/2, qp=1.6, r0=0.5)


#define our island grids...
Nκ = 500
Nᾱ = 100
Nφ = 10
#κgrid = LinRange(0, 8, Nκ);
κgrid = LinRange(0, 1.5, Nκ);
ᾱgrid = LinRange(-π, π, Nᾱ+1)[1:end-1];
φgrid = LinRange(0, 2π, Nφ+1)[1:end-1];

ϕ_isl  = create_island_mode(κgrid, ᾱgrid, φgrid, true);

surface(ᾱgrid, κgrid, real.(ϕ_isl[:, :, 1]))
contourf(ᾱgrid, κgrid, real.(ϕ_isl[:, :, 1]), color=:turbo)

hline!([0, 2π], [1, 1], color=:black, linewidth=5, legend=false)
#now we convert this island mode into regular coordinates.


isl_itp = interpolate((κgrid, ᾱgrid, φgrid), ϕ_isl, Gridded(Linear(Periodic())));

isl_ext = extrapolate(isl_itp, Periodic());

Nr = 500
Nθ = 100
Nζ = 10

ϕ_tor = zeros(ComplexF64, Nr, Nθ, Nζ);

rgrid = LinRange(0, 1, Nr);
θgrid = LinRange(-π, π, Nθ+1)[1:end-1];
ζgrid = LinRange(0, 2π, Nζ+1)[1:end-1];

κvals = zeros(Nr, Nθ, Nζ);
ᾱvals = zeros(Nr, Nθ, Nζ);
φvals = zeros(Nr, Nθ, Nζ);

#store the sepratrix location.
sep1 = zeros(Nθ)
sep2 = zeros(Nθ)

ψsep1 = zeros(Nθ)
ψsep2 = zeros(Nθ)


#this is going to assume ζ=0 for simplicity, should generalise though!
for i in 1:Nθ
    α = θgrid[i]
    #unsure why this needs to be isl.A + isl.A -> should be negative???
    res = sqrt(abs(-2 * isl.q0^2 / isl.qp * (isl.A + isl.A * cos(isl.m0 * α))))

    ψsep1[i] = -res + isl.r0
    ψsep2[i] = res + isl.r0
end




#function find_ind(val, grid)

#    return argmin(abs.(grid .- val))
#end

ψgrid = 1/2 .* rgrid .^ 2
#ψgrid = rgrid
for (i, ψ) in enumerate(ψgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)

    κ, ᾱ, φ = toroidal_to_isl(ψ, θ, ζ, isl)

    #κind = find_ind(κgrid, κ)
    #ᾱind = find_ind(ᾱgrid, mod(ᾱ, 2π))
    #φind = find_ind(φgrid, φ)

    #ϕ_tor[i, j, k] = ϕ_isl[κind, ᾱind, φind] 
    ϕ_tor[i, j, k] = isl_ext(sqrt(κ), ᾱ, φ)

    κvals[i, j, k] = κ
    #ᾱvals[i, j, k] = mod(ᾱ, 2π)
    ᾱvals[i, j, k] = ᾱ
    φvals[i, j, k] = φ
end

#outside case seems to be working ish, mode is weirdly spiky though??? unclear why??? -> perhaps because of resolution???


#this looks ok, as long as we mod ᾱ but the contour of ᾱ is still pretty cooked.
contourf(θgrid, rgrid, real.(ϕ_tor[:, :, 1]), levels=50, color=:turbo)

scatter!(θgrid, ψsep1, color=:black, markersize=2, legend=false)
scatter!(θgrid, ψsep2, color=:black, markersize=2, legend=false)

surface(θgrid, rgrid, real.(ϕ_tor[:, :, 1]), levels=20, color=:turbo)

#maybe the island should looks antisymmetric in radial coords?
#do we see that in our code?
#we may have to swap to flux........... fk me.

#using ψ not r fixes up down symmetry problemo, still have left right symmetry problemo.
#so I guess we need to work in flux... which isn't going to suit our larger purpose very well. as we will have an island defined in terms of r no ψ that we need to map...
#this seems to only work if I have a uniform ψ grid...

#not rotating? Should it be???
#these have got to be wrong. unsure what is going on though!
contourf(θgrid, rgrid, ᾱvals[:, :, 1], levels=50, color=:turbo)
contourf(θgrid, rgrid, κvals[:, :, 1], levels=50, color=:turbo)
contourf(θgrid, rgrid, φvals[:, :, 1], levels=50, color=:turbo)



#now, can we convert this back???

#define our island grids...
Nκ2 = 500
Nᾱ2 = 100
Nφ2 = 10
#κgrid2 = LinRange(0, 8, Nκ2);
κgrid2 = LinRange(0, 1.5, Nκ2);
ᾱgrid2 = LinRange(0, 2π, Nᾱ2+1)[1:end-1];
φgrid2 = LinRange(0, 2π, Nφ2+1)[1:end-1];

ϕ_isl2 = zeros(ComplexF64, Nκ2, Nᾱ2, Nφ2);

θvals = zeros(Nκ2, Nᾱ2, Nφ2);
rvals = zeros(Nκ2, Nᾱ2, Nφ2);
ζvals = zeros(Nκ2, Nᾱ2, Nφ2);


for (i, κ) in enumerate(κgrid2), (j, ᾱ) in enumerate(ᾱgrid2), (k, φ) in enumerate(φgrid2)

    r, θ, ζ = isl_to_toroidal(κ, ᾱ, φ, isl)

    θvals[i, j, k] = mod(θ, 2π)
    rvals[i, j, k] = r
    ζvals[i, j, k] = ζ

    rind = find_ind(rgrid, r)
    θind = find_ind(θgrid, mod(θ, 2π),)
    ζind = find_ind(ζgrid, ζ)

    ϕ_isl2[i, j, k] = ϕ_tor[rind, θind, ζind]
end

#sort of works, not super clear though...
contourf(ᾱgrid2, κgrid2, real.(ϕ_isl2[:, :, 1]), levels=50, color=:turbo)

hline!([0, 2π], [1, 1], color=:black, linewidth=5, legend=false)

contourf(ᾱgrid2, κgrid2, θvals[:, :, 1], levels=50, color=:turbo)
contourf(ᾱgrid2, κgrid2, rvals[:, :, 1], levels=50, color=:turbo)
contourf(ᾱgrid2, κgrid2, ζvals[:, :, 1], levels=50, color=:turbo)


#testing the fourier transform
using FFTW
using Printf

ϕ_isl_fft = fft(ϕ_isl, [2, 3])

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
    plot!(κgrid, real.(ϕ_isl_fft[:, i, j]), label=@sprintf("(%s, %s)", m, n))
end
#just the two zero mode. as expected...
#again just the (2, 0) mode.
#and again.
display(p)



ϕ_tor_fft = fft(ϕ_tor, [2, 3])

p = plot()
for i in 1:Nθ, j in 1:Nζ
    if i > Nθ / 2
        m = i - Nθ - 1
    else
        m = i-1
    end
    if j > Nζ/2
        n = j - Nζ - 1
    else
        n = j-1
    end

    plot!(rgrid, real.(ϕ_tor_fft[:, i, j]), label=@sprintf("(%s, %s)", m, n))
end
#inside global
#main mode is (0, 0)
#then (3, 2), (6, 4), (9, 6) etc.
#inside global
#main mode is (0, 0) still, but not as drastically.
#still see predominatly (3, 2) multiples
#but now there is a large influence of random nonsense.
#outside cont
#basically just giberish!
display(p)



ϕ_isl2_fft = fft(ϕ_isl2, [2, 3])

p = plot()
for i in 1:Nᾱ2, j in 1:Nφ2
    if i > Nᾱ2/2
        m = i - Nᾱ2 - 1
    else
        m = i-1
    end
    if j > Nφ2/2
        n = j - Nφ2 - 1
    else
        n = j-1
    end
    plot!(κgrid, real.(ϕ_isl2_fft[:, i, j]), label=@sprintf("(%s, %s)", m, n))
end
#(2, 0) and (-2, 0) are the dominant mode, but clearly not the only one
#also have (1, 0) and (3, 0)
#same for inside_cont case.
#outside cont
#still see dominat mode as expected, but getting a lot more gibbereish than before.
display(p)