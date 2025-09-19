#may need a new version that considers a (1, 1) island.
#as that is the most extrame example
#we may want to also try the (2, 1) island again, if (1, 1) is cooked af due to (1, 1) branch being a bit cooked.
using MID
using MIDViz
using MIDCantori
using Plots; plotlyjs()
using Plots; gr()
#%%
R0 = 1.0
#this seems like an adequate choice.
isl = init_island(m0=1, n0=-1, w=0.2, qp=1.0, ψ0=1/2, flux=true)
geo = init_geo(R0=R0)
prob = init_problem(q=island_q, met=:cylinder, geo=geo, isl=isl, type=:flux)
prob.isls[1]
M = 32
N = 8
#%%
ψgrid = init_grid(type=:rc, N=100)
θgrid = init_grid(type=:as, N=4, start=0)
ζgrid = init_grid(type=:as, N=4, start=-3)
grids = init_grids(ψgrid, θgrid, ζgrid)
ωa = MID.Solve.analytical_continuum(prob, grids)
continuum_plot(ωa)
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
rats1 = lowest_rationals(8, prob.q(0.8)[1], prob.q(1.0)[1])
#gl1 = surface_guess(rats1, prob.q)
gl1 = 0.9*ones(length(rats1))
surfs1 = construct_surfaces(rats1, gl1, prob, M=M, N=N);
plot_surfs(surfs1)
#%%
#now this is the spicy one!
rats2 = lowest_rationals(10, prob.q(0.6)[1], prob.q(0.8)[1])
gl2 = surface_guess(rats2, prob.q)
surfs2 = construct_surfaces(rats2, gl2, prob, M=M, N=N);
plot_surfs(surfs2)
#%%
rats3 = lowest_rationals(15, prob.q(0.4)[1], prob.q(0.6)[1])
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
#not needed for the (1, 1) island.
rats6 = [(29, 15), (5, 4), (31, 30), (51, 50)]
gl6 = surface_guess(rats6, prob.q)
surfs6 = construct_surfaces(rats6, gl6, prob, M=M, N=N);
plot_surfs(surfs6)
#%%
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5, surfs6);
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5);
plot_surfs(curr_surfs)
#%%
using JLD2
save_object("/Users/matt/phd/MID/data/surfaces/cantori_island/w2_surfs.jld2", curr_surfs)
