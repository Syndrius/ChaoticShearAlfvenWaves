#can we use the cantori set up for an island case?
using MID
using MIDViz
using MIDCantori
using Plots; plotlyjs()
using Plots; gr()
#%%

k10 = 0.0010 #0.103
k12 = 0.0012 #0.113
k17 = 0.0017 #0.1347
k20 = 0.0020 #w=0.146
k24 = 0.0024 #w=0.16
k30 = 0.0030 #w=0.1789
k36 = 0.0036 #w=0.1959

k = k36

R0 = 1.0
#this should match the cantori_q set up perfectly.
#noice, they look to match.
#isl = init_island(m0=3, n0=-2, w=0.1, qp=9/8, ψ0=2/3, flux=true)
isl = init_island(m0=3, n0=-2, A=k/3, qp=9/8, ψ0=2/3, flux=true)
geo = init_geo(R0=R0)
prob = init_problem(q=island_q, met=:cylinder, geo=geo, isl=isl, type=:flux)
prob.isls[1]
M = 32
N = 8
#%%
islc = init_island(m0=3, n0=-2, flux=true, A=0.0012/3)
probc = init_problem(q=cantori_q, met=:cylinder, geo=geo, isl=islc, type=:flux)
#%%

Ntraj = 150;
ψlist = collect(LinRange(0.005, 0.995, Ntraj));
#rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

poincare_plot(prob, Nlaps, ψlist, zeros(Ntraj))
#poincare_plot(probc, Nlaps, ψlist, zeros(Ntraj))

#%%
#now we just need to generate some surfaces and re-run the ol tests.
rats1 = lowest_rationals(10, prob.q(0.8)[1], prob.q(1.0)[1])
#gl1 = surface_guess(rats1, prob.q)
gl1 = 0.9*ones(length(rats1))
surfs1 = construct_surfaces(rats1, gl1, prob, M=M, N=N);
plot_surfs(surfs1)
#%%
#now this is the spicy one!
rats2 = lowest_rationals(11, prob.q(0.6)[1], prob.q(0.8)[1])
gl2 = surface_guess(rats2, prob.q)
surfs2 = construct_surfaces(rats2, gl2, prob, M=M, N=N);
plot_surfs(surfs2)
#%%
rats3 = lowest_rationals(13, prob.q(0.4)[1], prob.q(0.6)[1])
gl3 = surface_guess(rats3, prob.q)
surfs3 = construct_surfaces(rats3, gl3, prob, M=M, N=N);
plot_surfs(surfs3)
#%%
rats4 = lowest_rationals(15, prob.q(0.2)[1], prob.q(0.4)[1])
gl4 = surface_guess(rats4, prob.q)
surfs4 = construct_surfaces(rats4, gl4, prob, M=M, N=N);
plot_surfs(surfs4)
#%%
rats5 = lowest_rationals(19, prob.q(0.0)[1], prob.q(0.2)[1])
gl5 = surface_guess(rats5, prob.q)
surfs5 = construct_surfaces(rats5, gl5, prob, M=M, N=N);
plot_surfs(surfs5)
#%%
rats6 = [(29, 15), (5, 4), (31, 30), (51, 50)]
gl6 = surface_guess(rats6, prob.q)
surfs6 = construct_surfaces(rats6, gl6, prob, M=M, N=N);
plot_surfs(surfs6)
#%%
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5, surfs6);
plot_surfs(curr_surfs)
#%%
using JLD2
save_object("/Users/matt/phd/MID/data/surfaces/cantori_island/isl_k36_surfs.jld2", curr_surfs)
curr_surfs
