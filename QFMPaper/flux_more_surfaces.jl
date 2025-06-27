
using MID
using MIDViz
using JLD2
using Plots; plotlyjs()
using Plots; gr()
#%%
# may cause issues with the secondary islands!
#k = 0.001 could be a base case of practically no chaos!
k_min = 0.0012 #ideally we could get this theoretically!

k_mid = 0.0016

k_max = 0.0019 #above this the qfm surfaces inside the chaotic region get cooked, this may still be too much though!

k = k_min

geo = init_geo(R0=4.0)
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=5, n0=-3, A=k/5, flux=true)

isls = MID.Geometry.IslandT[isl1, isl2]

prob = init_problem(geo=geo, q=MID.Equilibrium.cyl_qfm_q, isls=isls, met=:cylinder, type=:flux)

M = 64
N = 32

#%%

Ntraj = 300;
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
ψlist = collect(LinRange(0.4, 0.8, Ntraj));
Nlaps = 500;

#would be good if we could plot the poincare with like a gradient of colours to see where they end up.
poincare_plot(prob, Nlaps, Ntraj, ψlist, ylimits=(0.5, 0.67))

#%%
rats1 = lowest_rationals(11, prob.q(0.0)[1], prob.q(1.0)[1])
gl1 = surface_guess(rats1, prob.q)
#changing these numbers doesn't really help remove the spikes
surfs1 = construct_surfaces(rats1, gl1, prob, M=M, N=N);
plot_surfs(surfs1)
#%%
rats2 = [(16, 15), (29, 15), (39, 20), (30, 29)]
#rats2 = [(30, 29), (23, 12), (39, 20)]
#rats2 = [(23, 12), (16, 15), (29, 15)]
#rats2 = [(39, 20)]
gl2 = surface_guess(rats2, prob.q)
surfs2 = construct_surfaces(rats2, gl2, prob, M=M, N=N);
plot_surfs(surfs2)
#%%
curr_surfs = vcat(surfs1, surfs2);
plot_surfs(curr_surfs)
#%%
save_object("/Users/matt/phd/MID/data/surfaces/flux_qfm_min_surfs.jld2", curr_surfs)
