

using MID
using JLD2
#dgl is now fine after fixing the bugs, the interpolated derivatives are slow to converge but they do!
#error is also minimal.
#dgu is essentially the same!

surfs = load_object("data/qfm/benchmark_surfs.jld2");
#plot_surfs(surfs)
surf_itp = MID.QFM.create_surf_itp(surfs);

function f(r, θ, ζ)

    #peak efficiency
    CT = MID.QFM.CoordTsfmT()
    tor_met = MID.Geometry.MetT()
    qfm_met = MID.Geometry.MetT()

    MID.QFM.coord_transform!(r, θ, ζ, CT, surf_itp)
    MID.Geometry.toroidal_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], 4.0)

    MID.QFM.met_transform!(tor_met, qfm_met, CT)

    return qfm_met.gu
    #return CT.JM
end
function df(r, θ, ζ)

    #peak efficiency
    CT = MID.QFM.CoordTsfmT()
    tor_met = MID.Geometry.MetT()
    qfm_met = MID.Geometry.MetT()

    MID.QFM.coord_transform!(r, θ, ζ, CT, surf_itp)
    MID.Geometry.toroidal_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], 4.0)

    MID.QFM.met_transform!(tor_met, qfm_met, CT)

    return qfm_met.dgu
end
function dfdr(r, θ, ζ, Δr)
    return @. (f(r+Δr, θ, ζ) - f(r - Δr, θ, ζ)) / (2 * Δr)
    return @. (-f(r+2Δr, θ, ζ) + 8 * f(r + Δr, θ, ζ) - 8 * f(r-Δr, θ, ζ) + f(r - 2*Δr, θ, ζ)) / (12 * Δr)
end
function dfdθ(r, θ, ζ, Δθ)
    return @. (f(r, θ+Δθ, ζ) - f(r, θ-Δθ, ζ)) / (2 * Δθ)
    return @. (-f(r, θ+2Δθ, ζ) + 8 * f(r, θ+Δθ, ζ) - 8 * f(r, θ-Δθ, ζ) + f(r, θ-2Δθ,  ζ)) / (12 * Δθ)
end
function dfdζ(r, θ, ζ, Δζ)
    return @. (f(r, θ, ζ+Δζ) - f(r, θ, ζ-Δζ)) / (2 * Δζ)
    return @. (-f(r, θ, ζ+2Δζ) + 8 * f(r, θ, ζ+Δζ) - 8 * f(r, θ, ζ-Δζ) + f(r, θ,  ζ-2Δζ)) / (12 * Δζ)
end
#%%

hlist = [0.1, 0.05, 0.03, 0.01, 0.005, 0.003, 0.001, 0.0005, 0.0003, 0.0001]
Nh = length(hlist)
Nr = 11
Nθ = 11
Nζ = 11
rvals = LinRange(0.5, 0.6, Nr)
θvals = LinRange(0.2, 6.0, Nθ)
ζvals = LinRange(0.2, 6.0, Nζ)

#each deriv is treated separatly!
error = zeros(Nh, Nr, Nθ, Nζ, 3, 3, 3);

for (i, h) in enumerate(hlist), (j, r) in enumerate(rvals), (k, θ) in enumerate(θvals), (l, ζ) in enumerate(ζvals)
    an = df(r, θ, ζ)
    #display((r, θ, ζ))

    error[i, j, k, l, 1, :, :] = abs.(an[:, :, 1] .- dfdr(r, θ, ζ, h))
    error[i, j, k, l, 2, :, :] = abs.(an[:, :, 2] .- dfdθ(r, θ, ζ, h))
    error[i, j, k, l, 3, :, :] = abs.(an[:, :, 3] .- dfdζ(r, θ, ζ, h))
end
#%% 

avg_error = zeros(Nh, 3, 3, 3);

for i in 1:Nh, j in 1:3, k in 1:3, l in 1:3
    #avg_error[i, j, k, l] = maximum(error[i, :, :, :, j, k, l])
    avg_error[i, j, k, l] = mean(error[i, :, :, :, j, k, l])
end
#%%
lh = log.(hlist);
le = log.(avg_error);

p = plot()#ylims=(0, 5))
for i in 1:3, j in 1:3, k in 1:3
    #plot!(-lh, le[:, i, j, k], label="d"*string(i)*"g_"*string(j) * string(k))
    plot!(-lh, avg_error[:, i, j, k], label="d"*string(i)*"g_"*string(j) * string(k))
end
display(p)


#%%
grad = zeros(Nh-1, 3, 3, 3);

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

function tor_test(r, θ, ζ)

    #peak efficiency
    CT = MID.QFM.CoordTsfmT()
    tor_met = MID.Geometry.MetT()
    qfm_met = MID.Geometry.MetT()

    MID.QFM.coord_transform!(r, θ, ζ, CT, surf_itp)
    MID.Geometry.toroidal_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], 4.0)

    MID.QFM.met_transform!(tor_met, qfm_met, CT)

    #display(met.gu .- t1)

    display(CT.coords)
    display(CT.JM[1, 3])

end

tor_test(0.1, 0.2, 0.3)


#%%

rtest = 0.5
θtest = 0.4
ζtest = 0.3
htest = 0.0001
dgldr = df(rtest, θtest, ζtest)[:, :, 1]
fdgldr = dfdr(rtest, θtest, ζtest, htest)
dgldθ = df(rtest, θtest, ζtest)[:, :, 2]
fdgldθ = dfdθ(rtest, θtest, ζtest, htest)
dgldζ = df(rtest, θtest, ζtest)[:, :, 3]
fdgldζ = dfdζ(rtest, θtest, ζtest, htest)
(abs.(dgldr .- fdgldr))
(abs.(dgldθ .- fdgldθ))
(abs.(dgldζ .- fdgldζ))

#%%

#testing the actual values of qfm_gl, for our un interp s case.
