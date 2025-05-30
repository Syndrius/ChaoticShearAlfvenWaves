#need to figure out how to map toroidal coords to qfm coords in order to do mapping.
#unfort I think we will need to do some interpolation which is suboptimal.
#there might be some clever fft way of doing this, but I am unsure


#Looks like it might be possible to find roots of the interpolated (r, θ, ζ) values to find the s such that r=r etc.

using MID
using MIDViz
using Plots; plotlyjs()
#%%


#%%

k1 = 0.02
geo = init_geo(R0=4.0)
isl1 = init_island(m0=3, n0=2, A=k1/3)

isls = [isl1]#, isl2, isl3]

prob = init_problem(geo=geo, q=low_shear_qfm_q, isls=isls)

#%%
#so (11, 7) does not work!, regardless of res
rats1 = lowest_rationals(10, prob.q(0.0)[1], prob.q(1.0)[1])
gl1 = surface_guess(rats1, prob.q)
#changing these numbers doesn't really help remove the spikes
@time surfs1 = construct_surfaces(rats1, gl1, prob, M=32, N=16);
plot_surfs(surfs1)
#%%
sgrid = init_grid(type=:rf, N=50, start=0.3, stop=0.9)
ϑgrid = init_grid(type=:af, N=6, pf=1)
φgrid = init_grid(type=:af, N=2, pf=-1)
#ϑgrid = init_grid(type=:as, N=2, start=1)
#φgrid = init_grid(type=:as, N=1, start=-1)
qfm_grids = init_grids(sgrid, ϑgrid, φgrid)

#%%
solver = init_solver(prob=prob, target=0.3, nev=100)
#%%
evals, ϕqfm, ϕqfmft = compute_spectrum_qfm(grids=qfm_grids, prob=prob, solver=solver, surfs=surfs1);
#%%
continuum_plot(evals)

qfm_ind = find_ind(evals, 0.287)
potential_plot(ϕqfmft, qfm_grids, qfm_ind, label_max=0.5)

#%%
rgrid = init_grid(type=:rf, N=50, start=0.3, stop=0.9);
θgrid = init_grid(type=:af, N=20);
ζgrid = init_grid(type=:af, N=5);
tor_grids = init_grids(rgrid, θgrid, ζgrid);
#%%
ϕ_tor = zeros(ComplexF64, rgrid.N, θgrid.N, ζgrid.N);
CT = MID.QFM.CoordTransformT();
surf_itp, sd = MID.QFM.create_surf_itp(surfs1);

rgrid, θgrid, ζgrid = MID.inst_grids(tor_grids);
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
MID.Mapping.map_qfm_to_tor(ϕ_tor, rgrid, θgrid, ζgrid, ϕqfmp, sgrid, ϑgp, φgp, CT, surf_itp, sd)
#%%
#actually looks perf I think, v similar mode, bit broader, and a bit more cooked, as expected!
#we are getting the symmetric modes i.e n=+- just as the island case, perhaps this is signifying the fft is a bit off?
#may need to think about that a wee bit.
#otherwise this is looking pretty good tbh!
#we probably should actually study the mapping from qfm->tor, as I would expect the mode structures to align with tor computations in limit.
ϕ_torft = fft(ϕ_tor, [2, 3]);
potential_plot(ϕ_torft, tor_grids, label_max=0.5)
#%%
#jokes this is cooked af lol, need higher res to see whats going on!
evals_tor, ϕtor, ϕtorft = compute_spectrum(grids=qfm_grids, prob=prob, solver=solver);

continuum_plot(evals_tor)
