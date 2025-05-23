#maybe this will look better with higher res in fff but by golly it looks shocking at the moment
using MID
using MIDViz
using JLD2
using Plots; plotlyjs()
#%%

rgrid = init_grid(type=:rc, N=100)
θgrid = init_grid(type=:as, N=4, start=1)
ζgrid = init_grid(type=:as, N=2, start=-2)
grids_cont = init_grids(rgrid, θgrid, ζgrid)
#%%

geo = init_geo(R0=10.0)
k = 0.002
k1 = k
k2 = k

isl = init_island(m0=2, n0=-1, A=k1/2)
#start with no islands
prob = init_problem(geo=geo, q=MID.Equilibrium.island_mode_21, met=:cylinder, isl=isl)
unprob = init_problem(geo=geo, q=MID.Equilibrium.island_mode_21, met=:cylinder)
#%%
evals_cont = compute_continuum(unprob, grids_cont);
#%%
continuum_plot(evals_cont, grids_cont)
#%%
surfs = load_object("/Users/matt/phd/MID/data/surfaces/island_surfs.jld2");
plot_surfs(surfs)
#%%
evals_cont_qfm = compute_continuum(prob, grids_cont, surfs);
#%%
continuum_plot(evals_cont_qfm, grids_cont)
#%%

rgrid = init_grid(type=:rf, N=30, start=0.1)
θgrid = init_grid(type=:af, N=10, pf=2)
ζgrid = init_grid(type=:as, N=5, start=-2)
grids = init_grids(rgrid, θgrid, ζgrid)
#%%
solver = init_solver(prob=prob, target=0.0)
#%%
evals_norm, ϕ_norm, ϕft_norm = compute_spectrum(prob=unprob, grids=grids, solver=solver);
evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=new_surfs);
evals_og, ϕ_og, ϕft_og = compute_spectrum(prob=prob, grids=grids, solver=solver);
#%%
continuum_plot(evals_norm)
continuum_plot(evals)
continuum_plot(evals_og)
#%%
function remove_surfs(rats, surfs)
    new_surfs = MID.QFM.QFMSurfaceT[]
    for surf in surfs
        if !(surf.rational in rats)
            push!(new_surfs, surf)
        end
    end
    return new_surfs
end

new_surfs = remove_surfs([(15, 7), (13, 6), (11, 5), (13, 7), (11, 6)], surfs);
