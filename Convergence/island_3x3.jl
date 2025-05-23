#jacobain is now working, magB is cooked, need to come back as there are probably earlier problems with B

using MID
#this file is for B/vectors rather than matrices
#this function will stand in for the different aspects we want to compute
#J and magB look good
function f(κ, ᾱ, ζ)

    k = 0.002
    #think this matches the island_mode_21, we will need further verification.
    isl = init_island(m0=2, n0=-1, A=k/2, r0=0.5, qp=2.0)
    isl = MID.Structures.inst_island(isl)
    #peak efficiency
    met = MID.Geometry.MetT()
    B = MID.Equilibrium.BFieldT()
    MID.Geometry.island_metric!(met, κ, ᾱ, ζ, 10.0, isl)

    MID.Equilibrium.compute_B_isl!(B, met, isl, κ, ᾱ, ζ)

    #isl1 = init_island(m0=3, n0=2, A=0.005)
    #isl2 = init_island(m0=1, n0=-4, A=0.0000)

    #MID.Equilibrium.compute_B!(B, met, qfm_benchmark_q, isl1, isl2, κ, ᾱ, ζ)

    
    return met.gl
end
function df(κ, ᾱ, ζ)

    k = 0.002
    #think this matches the island_mode_21, we will need further verification.
    isl = init_island(m0=2, n0=-1, A=k/2, r0=0.5, qp=2.0)
    isl = MID.Structures.inst_island(isl)
    #peak efficiency
    met = MID.Geometry.MetT()
    B = MID.Equilibrium.BFieldT()
    MID.Geometry.island_metric!(met, κ, ᾱ, ζ, 10.0, isl)

    MID.Equilibrium.compute_B_isl!(B, met, isl, κ, ᾱ, ζ)
    #isl1 = init_island(m0=3, n0=2, A=0.005)
    #isl2 = init_island(m0=1, n0=-4, A=0.0000)

    #MID.Equilibrium.compute_B!(B, met, qfm_benchmark_q, isl1, isl2, r, ᾱ, ζ)

    return met.dgl
end
function dfdκ(κ, ᾱ, ζ, Δκ)
    #return @. (f(κ+Δκ, ᾱ, ζ) - f(κ - Δκ, ᾱ, ζ)) / (2 * Δκ)
    return @. (-f(κ+2Δκ, ᾱ, ζ) + 8 * f(κ + Δκ, ᾱ, ζ) - 8 * f(κ-Δκ, ᾱ, ζ) + f(κ - 2*Δκ, ᾱ, ζ)) / (12 * Δκ)
end
function dfdᾱ(κ, ᾱ, ζ, Δᾱ)
    return @. (f(κ, ᾱ+Δᾱ, ζ) - f(κ, ᾱ-Δᾱ, ζ)) / (2 * Δᾱ)
end
function dfdζ(κ, ᾱ, ζ, Δζ)
    return @. (f(κ, ᾱ, ζ+Δζ) - f(κ, ᾱ, ζ-Δζ)) / (2 * Δζ)
end
#%%

hlist = [0.01, 0.005, 0.003, 0.001, 0.0005, 0.0003, 0.0001]
#hlist = [0.0003, 0.0001]
Nh = length(hlist)
Nκ = 21
Nᾱ = 21
Nζ = 21
κvals = LinRange(0.2, 0.8, Nκ)
ᾱvals = LinRange(0.2, 6.0, Nᾱ)
ζvals = LinRange(0.2, 6.0, Nζ)

#each deriv is treated separatly!
error = zeros(Nh, Nκ, Nᾱ, Nζ, 3, 3, 3);

for (i, h) in enumerate(hlist), (j, κ) in enumerate(κvals), (k, ᾱ) in enumerate(ᾱvals), (l, ζ) in enumerate(ζvals)
    an = df(κ, ᾱ, ζ)
    #display((κ, ᾱ, ζ))

    error[i, j, k, l, 1, :, :] = abs.(an[:, :, 1] - dfdκ(κ, ᾱ, ζ, h))
    error[i, j, k, l, 2, :, :] = abs.(an[:, :, 2] - dfdᾱ(κ, ᾱ, ζ, h))
    error[i, j, k, l, 3, :, :] = abs.(an[:, :, 3] - dfdζ(κ, ᾱ, ζ, h))
end
#% 

avg_error = zeros(Nh, 3, 3, 3);

for i in 1:Nh, j in 1:3, k in 1:3, l in 1:3
    #avg_error[i, j] = mean(error[i, :, :, :, j])
    avg_error[i, j, k, l] = maximum(error[i, :, :, :, j, k, l])
end
#%
lh = log.(hlist)
le = log.(avg_error)

p = plot()#ylims=(0, 5))
for i in 1:3, j in 1:3, k in 1:3
    #plot!(-lh, le[:, i], label="d"*string(i))
    #plot!(-lh, le[:, i, j, k], label="d"*string(i)*"g_"*string(j) * string(k))
    plot!(-lh, avg_error[:, i, j, k], label="d"*string(i)*"g_"*string(j) * string(k))
end
display(p)


#%%
#hopefully this is cooked, not the maths...
grad = zeros(Nh-1, 3, 3, 3)

for i in 1:Nh-1
    grad[i, :, :, :] = @. (le[i+1, :, :, :] - le[i, :, :, :]) / (lh[i+1] - lh[i])
end
p = plot()
for i in 1:3, j in 1:3, k in 1:3
    #plot!(-lh, le[:, i, j, k])
    plot!(-lh[1:end-1], grad[:, i, j, k], label="d"*string(i)*"g_"*string(j) * string(k))
end
display(p)



#%%
#yikes, something is wrong with dgu...

function tor_test(κ, ᾱ, ζ)

    #peak efficiency
    met = MID.Geometry.MetT()
    MID.Geometry.toroidal_metric!(met, r, ᾱ, ζ, 4.0)
    t1 = inv(met.gl)

    #display(met.gu .- t1)

    test = zeros(3, 3, 3)
    for i in 1:3
        test[:, :, i] = -t1 * met.dgl[:, :, i] * t1
    end
    #display(met.dgu .- test)

    display(met.gu)
    display(test[:, :, 1])
    display(met.dgu[:, :, 2])
    display(met.dgu[:, :, 3])
end

tor_test(0.1, 0.0, 0.0)


#%%

dgl = df(0.4, 0.4, 0.4)[:, :, 1]
fdgl = dfdr(0.4, 0.4, 0.4, 0.3)
log.(abs.(dgl .- fdgl))
