#testing our chaotic case in paralel, problems are showing on Gadi
using MID
using MIDParallel
using MIDViz
using Plots; plotlyjs()
#%%

geo = init_geo(R0=4.0)


k = 0.0003
isl1 = init_island(m0=6, n0=-5, A=1.0*k/6)
isl2 = init_island(m0=11, n0=-9, A=1.0*k/11)
isl3 = init_island(m0=9, n0=-7, A=1.0*k/6)
isl4 = init_island(m0=17, n0=-13, A=1.0*k/17)
isl5 = init_island(m0=13, n0=-11, A=1.0*k/13)
isl6 = init_island(m0=5, n0=-4, A=1.0*k/5)
isl7 = init_island(m0=13, n0=-10, A=1.0*k/13)
isl8 = init_island(m0=14, n0=-11, A=1.0*k/14)

isls = [isl1, isl2, isl3, isl4, isl5, isl6, isl7, isl8]

prob = init_problem(geo=geo, isls=isls, q=low_shear_q)
#%%

rgrid = init_grid(type=:rf, start=0.05, stop=0.95, N=40)
θgrid = init_grid(type=:as, start=0, N=6)
ζgrid = init_grid(type=:as, start=-4, N=4)
θgrid = init_grid(type=:af, pf=2, N=8)
ζgrid = init_grid(type=:af, pf=-2, N=4)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%
solver = init_solver(targets=[0.2, 0.3, 0.4], nev=200, prob=prob)
#%%

inputs_to_file(prob=prob, grids=grids, dir="data/", solver=solver)

#%%
par_post_process("data/")

#%%
evals = evals_from_file(dir="data/");

continuum_plot(evals)

ind = find_ind(evals, 0.35976)

ϕft = efunc_from_file(dir="data/", ind=ind);

potential_plot(ϕft, grids)
