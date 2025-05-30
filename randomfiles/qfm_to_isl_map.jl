#Doubt this will work in serial, but hopefully we can at least check the function runs without error
using MID
using MIDViz
using Plots; plotlyjs()
#%%

isl = init_island(m0=2, n0=-1, w=0.1, r0=0.5, qp=2.0)
geo = init_geo(R0=1000.0)
prob = init_problem(q=MID.Equilibrium.island_equiv_q, met=:cylinder, geo=geo, isl=isl)
#%%
rats1 = lowest_rationals(4, prob.q(0.0)[1], prob.q(0.8)[1])
#rats2 = lowest_rationals(10, prob.q(0.0)[1], prob.q(0.5)[1])

display(prob.q(0.0))
gl1 = surface_guess(rats1, prob.q)
#changing these numbers doesn't really help remove the spikes
@time surfs1 = construct_surfaces(rats1, gl1, prob, M=32, N=16);
plot_surfs(surfs1)
#%%
sgrid = init_grid(type=:rf, N=50, start=0.25, stop=0.75)
ϑgrid = init_grid(type=:af, N=6, pf=1)
φgrid = init_grid(type=:af, N=2, pf=-1)
#ϑgrid = init_grid(type=:as, N=2, start=1)
#φgrid = init_grid(type=:as, N=1, start=-1)
qfm_grids = init_grids(sgrid, ϑgrid, φgrid)

#%%
solver = init_solver(prob=prob, target=0.0, nev=100)
#%%
evals, ϕqfm, ϕqfmft = compute_spectrum_qfm(grids=qfm_grids, prob=prob, solver=solver, surfs=surfs1);
#%%
#3 options, none are particularly islandy, surprising right!
continuum_plot(evals)
qfm_ind = find_ind(evals, 0.08507)
potential_plot(ϕqfmft, qfm_grids, qfm_ind)
MIDViz.Plotting.contour_plot(ϕqfm, qfm_grids, ind=qfm_ind)
#%%
#nice thing about mapping qfm to island is that the boundaries in qfm being cooked does not matter!
κgrid = init_grid(type=:rf, N=50, start=0.0, stop=0.999);
ᾱgrid = init_grid(type=:af, N=20);
τgrid = init_grid(type=:af, N=5);
isl_grids = init_grids(κgrid, ᾱgrid, τgrid);
#%%
ϕ_isl = zeros(ComplexF64, κgrid.N, ᾱgrid.N, τgrid.N);
CT = MID.QFM.CoordTransformT();
surf_itp, sd = MID.QFM.create_surf_itp(surfs1);

κgrid, ᾱgrid, τgrid = MID.inst_grids(isl_grids);
sgrid, ϑgrid, φgrid = MID.inst_grids(qfm_grids);

ϕqfmp = zeros(ComplexF64, qfm_grids.x1.N, qfm_grids.x2.N+1, qfm_grids.x3.N+1);
ϕqfmp[:, 1:end-1, 1:end-1] = ϕqfm[qfm_ind, :, :, :];
ϕqfmp[:, :, end] = ϕqfmp[:, :, 1];
ϕqfmp[:, end, :] = ϕqfmp[:, 1, :];
ϕqfmp[:, end, end] = ϕqfmp[:, 1, 1];

ϑgp = LinRange(0, 2π, qfm_grids.x2.N+1);
φgp = LinRange(0, 2π, qfm_grids.x3.N+1);

#didn't error which is nice, this is pretty fkn slow, unsurprisingly.
#unsure what we can do about that tho!
#may just have to cop it
#perhaps run in parallel so we can at least do it somewhat quickly!
#I guess its not so bad, will just scale terribly, probably the same issues as the Hermite interpolation stuff
#ideally we would somehow map everything at once?
#no idea how though!
MID.Mapping.map_qfm_to_tor(ϕ_isl, κgrid, ᾱgrid, τgrid, ϕqfmp, sgrid, ϑgp, φgp, CT, surf_itp, sd)
#%%
#this doesn't really look like anything,
#hard to tell if this is a problem with the method or the starting mode.
#perhaps both
#we probably need higher res, and first map to toridal
#then to island to see whats what.
#none of the three really look anygood, but who knows tbh
ϕ_islft = fft(ϕ_isl, [2, 3]);
potential_plot(ϕ_islft, isl_grids, label_max=0.5)
