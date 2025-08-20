#hoepfully a more consistent appraoch to generating the surfaces
#here we will only consider the 21 case.
using MID
using MIDViz
using JLD2
using Plots; plotlyjs()
#%%
#cylindrical limit
geo = init_geo(R0=1.0)

#isl = init_island(m0=2, n0=-1, w=0.2, ψ0=0.5, qp=2.0, flux=true)
isl = init_island(m0=2, n0=-1, w=0.1, ψ0=0.5, qp=2.0, flux=true)

prob = init_problem(geo=geo, q=island_q, met=:cylinder, isl=isl, type=:flux)

prob.isls[1]
M = 32
N = 8
#%%
Ntraj = 150;
rlist = collect(LinRange(0.005, 0.995, Ntraj));
#rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

poincare_plot(prob, Nlaps, rlist, zeros(Ntraj))
#%%
rats1 = lowest_rationals(5, prob.q(0.8)[1], prob.q(1.0)[1])
#gl1 = surface_guess(rats1, prob.q)
gl1 = 0.9*ones(length(rats1))
surfs1 = construct_surfaces(rats1, gl1, prob, M=M, N=N);
plot_surfs(surfs1)
#%%
rats2 = lowest_rationals(7, prob.q(0.6)[1], prob.q(0.8)[1])
gl2 = surface_guess(rats2, prob.q)
surfs2 = construct_surfaces(rats2, gl2, prob, M=M, N=N);
plot_surfs(surfs2)
#%%
#this will be the spiciest one!
rats3 = lowest_rationals(7, prob.q(0.4)[1], prob.q(0.6)[1])
gl3 = surface_guess(rats3, prob.q)
surfs3 = construct_surfaces(rats3, gl3, prob, M=M, N=N);
plot_surfs(surfs3)
#%%
rats4 = lowest_rationals(10, prob.q(0.2)[1], prob.q(0.4)[1])
gl4 = surface_guess(rats4, prob.q)
surfs4 = construct_surfaces(rats4, gl4, prob, M=M, N=N);
plot_surfs(surfs4)
#%%
rats5 = lowest_rationals(15, prob.q(0.0)[1], prob.q(0.2)[1])
gl5 = surface_guess(rats5, prob.q)
surfs5 = construct_surfaces(rats5, gl5, prob, M=M, N=N);
plot_surfs(surfs5)
#%%
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5);
plot_surfs(curr_surfs)
#%%
#save_object("/Users/matt/phd/MID/data/surfaces/island/21a_surfs.jld2", curr_surfs)
save_object("/Users/matt/phd/MID/data/surfaces/island/w1_surfs.jld2", curr_surfs)
#save_object("/Users/matt/phd/MID/data/surfaces/island/w2_surfs.jld2", curr_surfs)
