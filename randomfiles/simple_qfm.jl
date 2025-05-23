#nothing in qfm is working properly, think we need to try a simpler case of either a cylinder or a slab.
#first we will test if a slab will work
#cant really remember how to deal with the slab normalisation stuff so I think we need will start with a cylindrical case.
#although that is still annoying as we have to deal with the axis then
#we may need to consider our orignal case in cylindrical limit.

using MID
using MIDViz
#%%

function simple_qfm_q(r::Float64)

    #so q varies form 1 to 2 linearly
    return 1 + r, 1
end

geo = init_geo(R0=100.0)
k = 0.005
k1 = k
k2 = k

isl1 = init_island(m0=4, n0=-3, A=k1/4)
isl2 = init_island(m0=3, n0=-2, A=k2/3)
isl3 = init_island(m0=5, n0=-3, A=k1/5)
isls = [isl1, isl2, isl3]
#start with no islands
prob = init_problem(geo=geo, q=simple_qfm_q, met=:cylinder, isls=isls)
unprob = init_problem(geo=geo, q=simple_qfm_q, met=:cylinder)
#%%

Ntraj = 150;
rlist = collect(LinRange(0.005, 0.995, Ntraj));
#rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

poincare_plot(prob, Nlaps, Ntraj, rlist)
#%%
rgrid_cont = init_grid(type=:rc, N = 100)
θgrid_cont = init_grid(type=:as, N=4, start=1)
ζgrid_cont = init_grid(type=:as, N=2, start=-2)

grids_cont = init_grids(rgrid_cont, θgrid_cont, ζgrid_cont)
#%%
evals_cont = compute_continuum(unprob, grids_cont);
#%%
continuum_plot(evals_cont, grids_cont)

#%%
rgrid = init_grid(type=:rf, N = 100)
θgrid = init_grid(type=:as, N=4, start=1)
ζgrid = init_grid(type=:as, N=2, start=-2)

grids = init_grids(rgrid, θgrid, ζgrid)

solver = init_solver(full_spectrum=true, prob=unprob)
#%%
evals, ϕ, ϕft = compute_spectrum(prob=unprob, grids=grids, solver=solver);
#%%

continuum_plot(evals)
#%%
#now create the qfm surfaces
#looks like there is to much chaos, surfaces in chaotic region cannot found
rats1 = lowest_rationals(11, prob.q(0.0)[1], prob.q(1.0)[1])
gl1 = surface_guess(rats1, prob.q)
#changing these numbers doesn't really help remove the spikes
surfs1 = construct_surfaces(rats1, gl1, prob, M=32, N=16);
