
#same as before but we will use Axel's transformation.
#naturally this will only work inside the island.

#so this is basically identical, but also gives some cooked stuff for sin waves, may be a hint at sol
#still seams like the transformation to toroidal is fine, just the reverse is not.

using Elliptic
using FFTW
using MID
using Plots; plotlyjs()
using Interpolations


function global_island_mode_dm(κ, βs, φ, isl)

    w = 0.1

    #κ = κ^2
    κ = sqrt(κ)
    if κ > 1 
        return 0
        #return 0, 0, 0, 0
    end

    K = Elliptic.K(κ) 
    sinβ = Elliptic.Jacobi.sn(4*K * βs / 2π, κ)
    cosβ = Elliptic.Jacobi.cn(4*K * βs / 2π, κ)

    r = sqrt(isl.ψ0*2  + w*κ * cosβ)

    α = asin(κ * sinβ)

    θ = (2 * α - isl.n0 * φ)/isl.m0

    ζ = φ


    α2 = 1/2 * (isl.m0*θ + isl.n0*ζ)

    

    r0 = sqrt(isl.ψ0*2)

    w = 0.1
    #srs confusion about this...
    #β2 = atan(w*sin(α2) / (r^2-r0^2))
    β2 = atan(w*sin(α) , (r^2-r0^2))

    #κ = sqrt(1/w^2 * (r^2-r0^2)^2 + sin(α)^2)

    κ2 = sqrt(1/w^2 * (r^2-r0^2)^2 + sin(α2)^2)^2

    if κ < 1
        βs2 =  π / (2 * Elliptic.K(κ2)) * Elliptic.F(β2, κ2)
    else
        #return 0
        return 0, 0, 0, 0
    end

    m = 2

    

    #return κ2, βs2, r, θ

    #return (-4*(κ-0.5)^2+1)*exp(1im * m * βs)
    #this perhaps shows the problemo...
    return (-4*(κ2-0.5)^2+1)*sin(m * βs2)
    #imag part of this is also cooke emundo.
    #problemos galor.
    #return exp(-(κ-0.7)^2/0.01) * exp(1im * m * βs)
end




function global_island_mode_tor(r, θ, ζ, isl)

    α = 1/2 * (isl.m0*θ + isl.n0*ζ)

    

    r0 = sqrt(isl.ψ0*2)

    w = 0.1
    #srs confusion about this...
    β = atan(w*sin(α) / (r^2-r0^2))
    #β = atan(w*sin(α) , (r^2-r0^2))

    κ = sqrt(1/w^2 * (r^2-r0^2)^2 + sin(α)^2)

    #κ = sqrt(1/w^2 * (r^2-r0^2)^2 + sin(α)^2)^2

    if κ < 1
        βs =  π / (2 * Elliptic.K(κ)) * Elliptic.F(β, κ)
    else
        βs = 0
    end

    m = 2

    if κ > 1 
        return 0
    end

    #return (-4*(κ-0.5)^2+1)*exp(1im * m * βs)
    #this perhaps shows the problemo...
    return (-4*(κ-0.5)^2+1)*sin(m * βs)
    #imag part of this is also cooke emundo.
    #problemos galor.
    #return exp(-(κ-0.7)^2/0.01) * exp(1im * m * βs)
end

#note axel's kappa is squared version of ours.
function global_island_mode(κ, βs, φ)

    m = 2

    if κ > 1 
        return 0
    end

    return (-4*(κ-0.5)^2+1)*exp(1im * m * βs)
    #this perhaps shows the problemo...
    #return (-4*(κ-0.5)^2+1)*sin(m * βs)
    #imag part of this is also cooke emundo.
    #problemos galor.
    #return exp(-(κ-0.7)^2/0.01) * exp(1im * m * βs)


end

Nκ = 200
Nβs = 80
Nφ = 30

κgrid = LinRange(0, 1.5, Nκ)
βsgrid = LinRange(-π, π, Nβs+1)[1:end-1]
φgrid = LinRange(0, 2π, Nφ+1)[1:end-1]

ϕ_isl = zeros(ComplexF64, Nκ, Nβs, Nφ);

for (i, κ) in enumerate(κgrid), (j, βs) in enumerate(βsgrid), (k, φ) in enumerate(φgrid)

    ϕ_isl[i, j, k] = global_island_mode(κ, βs, φ)
end

contourf(βsgrid, κgrid, real.(ϕ_isl[:, :, 1]), color=:turbo, levels=100)
contourf(βsgrid, κgrid, imag.(ϕ_isl[:, :, 1]), color=:turbo)

hline!([0, 2π], [1, 1], color=:black, linewidth=5, legend=false)

isl = ContIslandT(m0=3, n0=2, qp=1.6, q0=3/2, A = 1e-4, ψ0 = 0.125)

function isl_to_tor_axel(κ, βs, φ, isl)

    w = 0.1

    #κ = κ^2
    #κ = sqrt(κ)
    if κ > 1 
        return 0, 0, 0
    end

    K = Elliptic.K(κ) 
    sinβ = Elliptic.Jacobi.sn(4*K * βs / 2π, κ)
    cosβ = Elliptic.Jacobi.cn(4*K * βs / 2π, κ)

    r = sqrt(isl.ψ0*2  + w*κ * cosβ)

    α = asin(κ * sinβ)

    θ = (2 * α - isl.n0 * φ)/isl.m0

    return r, θ, φ
end

function tor_to_isl_axel(r, θ, ζ, isl)

    α = 1/2 * (isl.m0*θ + isl.n0*ζ)

    

    r0 = sqrt(isl.ψ0*2)

    w = 0.1
    #srs confusion about this...
    #β = atan(w*sin(α) / (r^2-r0^2))
    β = atan(w*sin(α) , (r^2-r0^2))

    κ = sqrt(1/w^2 * (r^2-r0^2)^2 + sin(α)^2)

    #κ = sqrt(1/w^2 * (r^2-r0^2)^2 + sin(α)^2)^2

    if κ < 1
        βs =  π / (2 * Elliptic.K(κ)) * Elliptic.F(β, κ)
    else
        βs = 0
    end

    return κ, βs, ζ
end



#should actually be $r$ now.
Nr = 300
Nθ = 100
Nζ = 20

ϕ_tor = zeros(ComplexF64, Nr, Nθ, Nζ);
ϕ_tor_an = zeros(ComplexF64, Nr, Nθ, Nζ);

isl_itp = interpolate((κgrid, βsgrid, φgrid), ϕ_isl, Gridded(Linear(Periodic())));

isl_ext = extrapolate(isl_itp, Periodic());


rgrid = LinRange(0.0, 1, Nr);
θgrid = LinRange(-π, π, Nθ+1)[1:end-1];
ζgrid = LinRange(0, 2π, Nζ+1)[1:end-1];

for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)

    κ, βs, φ = tor_to_isl_axel(r, θ, ζ, isl)

    #κind = find_ind(κgrid, κ)
    #βsind = find_ind(βsgrid, rem(βs, π))
    #φind = find_ind(φgrid, rem(φ, π))

    #ϕ_tor[i, j, k] = #ϕ_isl[κind, βsind, φind] 
    #isl_itp(κ, βs, φ)
    ϕ_tor[i, j, k] = isl_ext(κ, βs, φ)
    ϕ_tor_an[i, j, k] = global_island_mode_tor(r, θ, ζ, isl)

    #κvals[i, j, k] = κ
    #ᾱvals[i, j, k] = mod(ᾱ, 2π)
    #ᾱvals[i, j, k] = ᾱ
    #φvals[i, j, k] = φ
end
#looks similar, but now we have the asymetry, as expected.
#this shows that our crummy interpolation method could be an actual srs problemo
contourf(θgrid, rgrid, real.(ϕ_tor[:, :, 11]), levels=100, color=:turbo)
contourf(θgrid, rgrid, real.(ϕ_tor_an[:, :, 11]), levels=100, color=:turbo)
contourf(θgrid, rgrid, imag.(ϕ_tor[:, :, 1]), levels=100, color=:turbo)



Nκ2 = 300
Nβs2 = 100
Nφ2 = 30

κgrid2 = LinRange(0, 1.5, Nκ2)
βsgrid2 = LinRange(-π, π, Nβs2+1)[1:end-1]
φgrid2 = LinRange(0, 2π, Nφ2+1)[1:end-1]


#tor_itp = LinearInterpolation((rgrid, θgrid, ζgrid), real.(ϕ_tor))#, BSpline(Linear()));
tor_itp = interpolate((collect(rgrid), collect(θgrid), collect(ζgrid)), ϕ_tor, Gridded(Linear(Periodic())))#, BSpline(Linear(Periodic())));

tor_ext = extrapolate(tor_itp, Periodic())#, Flat());


#tor_ext(25, 0.1, 0.1)
ϕ_isl2 = zeros(ComplexF64, Nκ2, Nβs2, Nφ2);
ϕ_isl2_an = zeros(ComplexF64, Nκ2, Nβs2, Nφ2);
ϕ_isl2_dm = zeros(ComplexF64, Nκ2, Nβs2, Nφ2);


for (i, κ) in enumerate(κgrid2), (j, βs) in enumerate(βsgrid2), (k, φ) in enumerate(φgrid2)

    r, θ, ζ = isl_to_tor_axel(κ, βs, φ, isl)

    rind = rind = find_ind(rgrid, r)
    θind = find_ind(θgrid, rem(θ, π))
    ζind = find_ind(ζgrid, rem(ζ, π))

    #ϕ_isl2[i, j, k] = ϕ_tor[rind, θind, ζind]
    ϕ_isl2[i, j, k] = tor_ext(r, θ, ζ)
    ϕ_isl2_an[i, j, k] = ϕ_tor_an[rind, θind, ζind]
    #ϕ_isl2_an[i, j, k] = ϕ_tor_an[rind, θind, ζind]
    ϕ_isl2_dm[i, j, k] = global_island_mode_dm(κ, βs, φ, isl)
end
#analytical one is significantly smoother, but it is still cooked af for peaks and troughs etc.
#may need to introduce some simple interpolation for the coordinate mapping...
#maybe we should come up with some arbitrary af function in regular coordinates and map that to island? See if we get something clear??? Like a liner function of kappa or something?

#this now seems to show that this is working well enough I guess.
contourf(βsgrid2, κgrid2, real.(ϕ_isl2[:, :, 1]), color=:turbo, levels=100)
contourf(βsgrid2, κgrid2, real.(ϕ_isl2_an[:, :, 1]), color=:turbo, levels=100)
contourf(βsgrid2, κgrid2, real.(ϕ_isl2_dm[:, :, 1]), color=:turbo, levels=100)
contourf(βsgrid2, κgrid2, imag.(ϕ_isl2[:, :, 1]), color=:turbo, levels=100)
contourf(βsgrid2, κgrid2, abs.(ϕ_isl2[:, :, 1]), color=:turbo, levels=100)
contourf(βsgrid, κgrid, real.(ϕ_isl[:, :, 1]), color=:turbo, levels=100)





using FFTW
using Printf


ϕ_isl_fft = fft(ϕ_isl, [2, 3]);
ϕ_tor_fft = fft(ϕ_tor, [2, 3]);
ϕ_isl2_fft = fft(ϕ_isl2, [2, 3]);



#think this looks stupid af because of the extended domain...
p = plot()
for i in 1:Nβs, j in 1:Nφ
    if i > Nβs/2
        m = i - Nβs - 1
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


#Ok, I think this is working to an adequate standard now, interesting that we get the lil (0, 0) tail at then end...
p = plot()
for i in 1:Nβs2, j in 1:Nφ2
    if i > Nβs2/2
        m = i - Nβs2 - 1
    else
        m = i-1
    end
    if j > Nφ2/2
        n = j - Nφ2 - 1
    else
        n = j-1
    end
    plot!(κgrid2, real.(ϕ_isl2_fft[:, i, j]), label=@sprintf("(%s, %s)", m, n))
end
#just the two zero mode. as expected...
#again just the (2, 0) mode.
#and again.
display(p)




#I think this tells us there is something wrong with the grids we are sampling from or the way we do mod
#I maybe stuff is extended beyond 2π?
contourf(βsgrid, κgrid, real.(ϕ_isl[:, :, 1]), color=:turbo, levels=100)
contourf(βsgrid2, κgrid2, real.(ϕ_isl2_dm[:, :, 1]), color=:turbo, levels=100)


klist = LinRange(0, 0.99999, 100)

#well fk me, julia's is built of k^2...
#this didn't seem to change any results tbh. but has to be a problemo
plot(klist, Elliptic.K.(klist.^2))
plot(klist, Elliptic.K.(klist.^2))

#think we should try a small κ m=0 case, and see what it looks like analytically.

Elliptic.K(0.5)



κ2vals = zeros(Nκ2, Nβs2, Nφ2);
βs2vals = zeros(Nκ2, Nβs2, Nφ2);
κvals = zeros(Nκ2, Nβs2, Nφ2);
βsvals = zeros(Nκ2, Nβs2, Nφ2);
rvals = zeros(Nκ2, Nβs2, Nφ2);
θvals = zeros(Nκ2, Nβs2, Nφ2);
for (i, κ) in enumerate(κgrid2), (j, βs) in enumerate(βsgrid2), (k, φ) in enumerate(φgrid2)

    κ2vals[i, j, k], βs2vals[i, j, k], rvals[i, j, k], θvals[i, j, k] = global_island_mode_dm(κ, βs, φ, isl)

    κvals[i, j, k] = κ
    βsvals[i, j, k] = βs

end

contourf(βsgrid2, κgrid2, real.(κ2vals[:, :, 1]), color=:turbo, levels=100)
contourf(βsgrid2, κgrid2, real.(βs2vals[:, :, 1]), color=:turbo, levels=100)
contourf(βsgrid2, κgrid2, real.(rvals[:, :, 1]), color=:turbo, levels=100)
contourf(βsgrid2, κgrid2, real.(θvals[:, :, 1]), color=:turbo, levels=100)
contourf(βsgrid2, κgrid2, real.(κvals[:, :, 1]), color=:turbo, levels=100)
contourf(βsgrid2, κgrid2, real.(βsvals[:, :, 1]), color=:turbo, levels=100)

