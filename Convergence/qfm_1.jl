

using MID
using JLD2
using Statistics
using Plots; plotlyjs()
#Jacobian seems to be fine now, did require changing our method of comptuation, unfort this means we are using the sqrt again,
#which may cause problemos later.
#dmag_B is also good.
#surfs = load_object("data/qfm/benchmark_surfs.jld2");
surfs = load_object("../data/qfm/benchmark_surfs.jld2");
#plot_surfs(surfs)
surf_itp = MID.QFM.create_surf_itp(surfs);
#this file is for B/vectors rather than matrices
#this function will stand in for the different aspects we want to compute
#J and magB look good
function f(r, θ, ζ)

    #peak efficiency
    CT = MID.QFM.CoordTsfmT()
    tor_met = MID.Geometry.MetT()
    qfm_met = MID.Geometry.MetT()
    tor_B = MID.Equilibrium.BFieldT()
    qfm_B = MID.Equilibrium.BFieldT()
    isl1 = init_island(m0=3, n0=2, A=0.005)
    isl2 = init_island(m0=1, n0=-4, A=0.0000)

    MID.QFM.coord_transform!(r, θ, ζ, CT, surf_itp)
    MID.Geometry.toroidal_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], 4.0)
    MID.Equilibrium.compute_B!(tor_B, tor_met, qfm_benchmark_q, isl1, isl2, CT.coords[1], CT.coords[2], CT.coords[3])

    MID.QFM.met_transform!(tor_met, qfm_met, CT)
    MID.QFM.B_transform!(tor_B, qfm_B, qfm_met, CT)

    return qfm_B.mag_B[1]
end
function df(r, θ, ζ)

    #peak efficiency
    CT = MID.QFM.CoordTsfmT()
    tor_met = MID.Geometry.MetT()
    qfm_met = MID.Geometry.MetT()
    tor_B = MID.Equilibrium.BFieldT()
    qfm_B = MID.Equilibrium.BFieldT()
    isl1 = init_island(m0=3, n0=2, A=0.005)
    isl2 = init_island(m0=1, n0=-4, A=0.0000)

    MID.QFM.coord_transform!(r, θ, ζ, CT, surf_itp)
    MID.Geometry.toroidal_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], 4.0)
    MID.Equilibrium.compute_B!(tor_B, tor_met, qfm_benchmark_q, isl1, isl2, CT.coords[1], CT.coords[2], CT.coords[3])

    MID.QFM.met_transform!(tor_met, qfm_met, CT)
    MID.QFM.B_transform!(tor_B, qfm_B, qfm_met, CT)

    return qfm_B.dmag_B
end
function dfdr(r, θ, ζ, Δr)
    return @. (f(r+Δr, θ, ζ) - f(r - Δr, θ, ζ)) / (2 * Δr)
    #return @. (-f(r+2Δr, θ, ζ) + 8 * f(r + Δr, θ, ζ) - 8 * f(r-Δr, θ, ζ) + f(r - 2*Δr, θ, ζ)) / (12 * Δr)
end
function dfdθ(r, θ, ζ, Δθ)
    return @. (f(r, θ+Δθ, ζ) - f(r, θ-Δθ, ζ)) / (2 * Δθ)
end
function dfdζ(r, θ, ζ, Δζ)
    return @. (f(r, θ, ζ+Δζ) - f(r, θ, ζ-Δζ)) / (2 * Δζ)
end
#%%

hlist = [0.01, 0.005, 0.003, 0.001, 0.0005, 0.0003, 0.0001]
hlist = [0.1, 0.05, 0.03, 0.01, 0.005, 0.003, 0.001, 0.0005, 0.0003, 0.0001]
Nh = length(hlist)
Nr = 11
Nθ = 11
Nζ = 11
rvals = LinRange(0.5, 0.6, Nr)
θvals = LinRange(0.2, 6.0, Nθ)
ζvals = LinRange(0.2, 6.0, Nζ)

#each deriv is treated separatly!
error = zeros(Nh, Nr, Nθ, Nζ, 3);

for (i, h) in enumerate(hlist), (j, r) in enumerate(rvals), (k, θ) in enumerate(θvals), (l, ζ) in enumerate(ζvals)
    an = df(r, θ, ζ)
    #display((r, θ, ζ))

    error[i, j, k, l, 1] = abs(an[1] - dfdr(r, θ, ζ, h))
    error[i, j, k, l, 2] = abs(an[2] - dfdθ(r, θ, ζ, h))
    error[i, j, k, l, 3] = abs(an[3] - dfdζ(r, θ, ζ, h))
end
#%% 

avg_error = zeros(Nh, 3);

for i in 1:Nh, j in 1:3
    avg_error[i, j] = mean(error[i, :, :, :, j])
    #avg_error[i, j] = maximum(error[i, :, :, :, j])
end
#%%
lh = log.(hlist)
le = log.(avg_error)

p = plot()#ylims=(0, 5))
for i in 1:3
    #plot!(-lh, le[:, i], label="d"*string(i))
    plot!(-lh, avg_error[:, i], label="d"*string(i))
end
display(p)


#%%
grad = zeros(Nh-1, 3)

for i in 1:Nh-1
    grad[i, :] = @. (le[i+1, :] - le[i, :]) / (lh[i+1] - lh[i])
end
p = plot()
for i in 1:3
    #plot!(-lh, le[:, i, j, k])
    plot!(-lh[1:end-1], grad[:, i], label="d"*string(i))
end
display(p)



#%%
#yikes, something is wrong with dgu...

function tor_test(r, θ, ζ)

    #peak efficiency
    met = MID.Geometry.MetT()
    MID.Geometry.toroidal_metric!(met, r, θ, ζ, 4.0)
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
