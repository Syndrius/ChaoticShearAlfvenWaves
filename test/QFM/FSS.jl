
sgrid = init_grid(type=:rf, N=30, start=0.4, stop=0.9)
ϑgrid = init_grid(type=:as, N=2, start=1)
φgrid = init_grid(type=:as, N=1, start=-1)
grids = init_grids(sgrid, ϑgrid, φgrid)

#%%

#solver = init_solver(nev=100, target = 0.3, prob=prob)
solver = init_solver(full_spectrum=true, prob=prob)

#%%

evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=surfs);

#continuum_plot(evals)

tae_ind = find_ind(evals, 0.2783)
tae_freq = evals.ω[tae_ind]

#potential_plot(ϕft, grids, tae_ind)

@test tae_ind == 58
@test tae_freq ≈ 0.2783 atol=0.01

