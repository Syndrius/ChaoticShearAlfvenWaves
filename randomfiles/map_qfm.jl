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
rats1 = lowest_rationals(10, prob.q(0.0)[1], prob.q(1.0)[1])
gl1 = surface_guess(rats1, prob.q)
@time surfs1 = construct_surfaces(rats1, gl1, prob, M=16, N=8);
plot_surfs(surfs1)
#%%
sgrid = init_grid(type=:rf, N=50, start=0.3, stop=0.9)
ϑgrid = init_grid(type=:af, N=6, pf=1)
φgrid = init_grid(type=:af, N=2, pf=-1)
#sgrid = init_grid(type=:rf, N=10, start=0.3, stop=0.9)
#ϑgrid = init_grid(type=:af, N=2, pf=1)
#φgrid = init_grid(type=:af, N=1, pf=-1)
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
rgrid = init_grid(type=:rf, N=100, start=0.35, stop=0.85);
θgrid = init_grid(type=:af, N=50);
ζgrid = init_grid(type=:af, N=10);
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
#think if we want to persist with this, we probably need to do it in parallel with like 100 cores.
#given each efunc takes closer to hours than minutes for practical cases.
#@time MID.Mapping.map_qfm_to_tor(ϕ_tor, rgrid, θgrid, ζgrid, ϕqfmp, sgrid, ϑgp, φgp, CT, surf_itp, sd)
#%%
#testing creating mapping structures once, so that we can very quickly map all of the efuncs.
#we will need the deriv for hermite stuff, but that can wait
coord_map = Array{Vector{Float64}}(undef, tor_grids.x1.N, tor_grids.x2.N, tor_grids.x3.N);


for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)
    #this is unfort a bit fked.
    #assuming this work ok though, this should be fine.
    coord_map[i, j, k] = MID.Mapping.tor_coords_to_qfm(r, θ, ζ, CT, surf_itp, sd)
end
#%%
@time MID.Mapping.qfm_efunc_to_tor(ϕ_tor, rgrid, θgrid, ζgrid, ϕqfmp, sgrid, ϑgp, φgp, coord_map)
#%%
#think we are running into problems due to the limited grid for qfm
#also seems to be a problmo if the qfm grid goes from 0.05 to 0.95, our tor grid should go from 0.08 to 0.92 or so
#pretty much won't be an issue for island case.
@time MID.Mapping.qfm_efunc_to_tor!(ϕ_tor, rgrid, θgrid, ζgrid, ϕqfm[qfm_ind, :, :, :, :], sgrid, ϑgrid, φgrid, coord_map)
#now lets see if we can make a Hermite map, given that the hermite things are computed at the same point every time
#this will be quite a large data structure unfort
#this is way to complicated for the effort, we would need to know each grid_ind and basis_ind.
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
#%%

display(coord_map[coord_map .< 0.0])
#now a comparison with the alternative way of mapping
#here we do not start with the mapped grids
#rather we create a new grid in the original coords that is much denser, hermite interpolate the eigenfunction
#then directly map that, should be much faster, but will result in a weird uneven grid.
#we may have to deal with that when we get to it.
#our grid instance is probably going to a problemo.
smgrid = init_grid(type=:rf, N=100, start=0.3, stop=0.9)
ϑmgrid = init_grid(type=:af, N=50, pf=1)
φmgrid = init_grid(type=:af, N=20, pf=-1)
#ϑgrid = init_grid(type=:as, N=2, start=1)
#φgrid = init_grid(type=:as, N=1, start=-1)
qfm_mgrids = init_grids(smgrid, ϑmgrid, φmgrid)


