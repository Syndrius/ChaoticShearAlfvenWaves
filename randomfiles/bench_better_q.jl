#this should not be called better_q, this is the new qfm q, we are generating the benchmark surfaces now.

using MID
using MIDViz
using JLD2

#%%
function q_prof(r::Float64)
    #chosen so (4, 3)=0.4, (3, 2)=0.6
    #sillier values are same but with flux surfaces shifted from 0, 1 to (0.05, 0.95)
    a = 6/5
    a = 539/450
    b = 5/6
    b = 8/9
    return a + b*r^2, 2 * b * r
end
#%%

rgrid = init_grid(type=:rf, N = 150, start=0.05, stop=0.95)
θgrid = init_grid(type=:as, start=-1, N = 6)
ζgrid = init_grid(type=:as, start=-5, N = 5)

grids = init_grids(rgrid, θgrid, ζgrid)

#%%

geo = init_geo(R0=4.0)
k = 0.001
isl = init_island(m0=4, n0=3, A = k/4)
prob = init_problem(geo=geo, q=q_prof, isl=isl)
unprob = init_problem(geo=geo, q=q_prof)

#%%

solver = init_solver(nev=200, targets=[0.23, 0.3, 0.37], prob=unprob)

#%%

evals_norm, ϕ_norm, ϕft_norm = compute_spectrum(prob=unprob, grids=grids, solver=solver);
#%%

continuum_plot(evals_norm)
ind_norm = find_ind(evals_norm, 0.2944)

potential_plot(ϕft_norm, grids, ind_norm)

#%%

Ntraj = 100;
#rlist = collect(LinRange(0.45, 0.65, Ntraj));
rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

#poincare_plot(qfm_benchmark_q, slab_to_plot, Nlaps, Ntraj, 0, 0.0, 0.0, 0.0, R0, isl, MID.Structures.no_isl, rlist)
poincare_plot(prob, Nlaps, Ntraj, rlist,  prob.geo.R0)#, filename=save_dir * "original_poincare.png");




chaos_surfs = load_object("/Users/matt/phd/MID/data/surfaces/chaos_surfs.jld2")
plot_surfs(chaos_surfs)
