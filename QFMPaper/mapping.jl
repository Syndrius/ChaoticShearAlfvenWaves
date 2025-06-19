
#testing how it looks when we map our solution from qfm space to cylindrical space
#hopefully the SAW will follow the cantori.
using MID
using MIDViz
using JLD2
using Plots; plotlyjs()
#%%

k1 = 0.00055
#k2 = 0.0003
#k3 = 0.0000
geo = init_geo(R0=4.0)
isl1 = init_island(m0=3, n0=-2, A=k1, flux=true)
isl2 = init_island(m0=5, n0=-3, A=k1, flux=true)
#isl3 = init_island(m0=8, n0=-5, A=k3/8) #unsure if we will want this one as well

isls = MID.Geometry.IslandT[isl1, isl2]

prob = init_problem(geo=geo, q=MID.Equilibrium.cyl_qfm_q, isls=isls, met=:cylinder, type=:flux)
#%%
curr_surfs = load_object("/Users/matt/phd/MID/data/cyl_qfm_surfaces/flux_surfaces.jld2");
#%%
sgrid = init_grid(type=:rf, N = 50, start=0.05, stop=0.95, sep1=0.5, sep2=0.66, frac=0.5)
ϑgrid = init_grid(type=:af, N = 6, pf=2) 
φgrid = init_grid(type=:af, N = 1, pf=-1)
grids = init_grids(sgrid, ϑgrid, φgrid)
#%%
solver = init_solver(nev=150, targets=[0.3], prob=prob)
#%%
evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=curr_surfs, deriv=true);
#%%

continuum_plot(evals, legend=false)#, n=-2)

ind = find_ind(evals, 0.2556594)

potential_plot(ϕft, grids, ind)

contour_plot(ϕ, grids, ind=ind)
#%%
rgrid = init_grid(type=:rf, N=100, start=0.15, stop=0.85)
θgrid = init_grid(type=:af, N=10)
ζgrid = init_grid(type=:af, N=5)

tor_grids = init_grids(rgrid, θgrid, ζgrid)

r_inst, θ_inst, ζ_inst = MID.inst_grids(tor_grids)
#%%
surf_itp, sd = MID.create_surf_itp(curr_surfs);
CT = MID.CoordTransformT();
coord_map = MID.Mapping.qfm_to_tor_coord_map(r_inst, θ_inst, ζ_inst, CT, surf_itp, sd);
#%%
ϕ_tor = zeros(ComplexF64, rgrid.N, θgrid.N, ζgrid.N);
s_inst, ϑ_inst, φ_inst = MID.inst_grids(grids);
#seems od to just pass the array shape in...
MID.Mapping.efunc_map!(ϕ_tor, rgrid.N, θgrid.N, ζgrid.N, ϕ[ind, :, :, :, :], s_inst, ϑ_inst, φ_inst, coord_map)

potential_plot(ϕ_tor, tor_grids)
#resolution is kind of garbage
#but this is pretty fkn sick
#defo shows what we want
#the interesting part will be to see if non-qfm approach can get this kind of shape to form.
contour_plot(ϕ_tor, tor_grids)
MIDViz.Plotting.overlay_surfs(curr_surfs)

