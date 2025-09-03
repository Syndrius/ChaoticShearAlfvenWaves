#note that we are not working with radius anymore!
#need to do a k1= 0.0011 case now, as that is a similar distance from the two criticals
#they are at ~0.0012, and 0.0014.
#the k2=0.0013 case is quite a good example now.
#may also want to change k3 = 0.0015, then we may not need k4.

using MID
using MIDViz
using JLD2
using Plots; plotlyjs()
using Plots; gr()
#%%

#our two surfaces have kc=0.0012 and kc = 0.0014.
#we have a case already with k = 0.0005, almost as a benchark, that is good,
#we also have a k = 0.0013 labelled as k2, which gives good results.
#we want k = 0.0015 as the maximum probably, although we may want to go higher if needed.
k11 = 0.0011 #this is quite the onset of chaos case.
#k13 = 0.0013 #already created the surfs for this.
k15 = 0.0015
k05 = 0.0005
k13 = 0.0013
k17 = 0.0017

k = k17
geo = init_geo(R0=1.0)
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)

isls = MID.Geometry.IslandT[isl1, isl2]

prob = init_problem(geo=geo, q=cantori_q, isls=isls, met=:cylinder, type=:flux)

M = 32
N = 8

#%%

Ntraj = 300;
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
ψlist = collect(LinRange(0.4, 0.8, Ntraj));
Nlaps = 500;

#would be good if we could plot the poincare with like a gradient of colours to see where they end up.
poincare_plot(prob, Nlaps, ψlist, zeros(Ntraj), ylimits=(0.5, 0.67))

#%%
rats1 = lowest_rationals(11, prob.q(0.0)[1], prob.q(1.0)[1])
gl1 = surface_guess(rats1, prob.q)
#changing these numbers doesn't really help remove the spikes
surfs1 = construct_surfaces(rats1, gl1, prob, M=M, N=N);
plot_surfs(surfs1)
#%%
rats2 = lowest_rationals(21, prob.q(0.96)[1], prob.q(1.0)[1])
#rats2 = [(16, 15), (29, 15), (39, 20), (30, 29)]
#rats2 = [(30, 29), (23, 12), (39, 20)]
#rats2 = [(23, 12), (16, 15), (29, 15)]
#rats2 = [(39, 20)]
gl2 = surface_guess(rats2, prob.q)
surfs2 = construct_surfaces(rats2, gl2, prob, M=M, N=N);
plot_surfs(surfs2)
#%%
rats3 = lowest_rationals(25, prob.q(0.0)[1], prob.q(0.15)[1])
gl3 = surface_guess(rats3, prob.q)
#changing these numbers doesn't really help remove the spikes
surfs3 = construct_surfaces(rats3, gl3, prob, M=M, N=N);
plot_surfs(surfs3)
#%%
rats4 = [(52, 51), (44, 43), (31, 30)]
gl4 = surface_guess(rats4, prob.q)
#changing these numbers doesn't really help remove the spikes
surfs4 = construct_surfaces(rats4, gl4, prob, M=M, N=N);
plot_surfs(surfs4)
#%%
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4);
plot_surfs(curr_surfs)
#%%
save_object("/Users/matt/phd/MID/data/surfaces/qfm/k17_surfs.jld2", curr_surfs)
