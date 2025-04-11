

rgrid = init_grid(type=:rf, N=15, start=0.4);
θgrid = init_grid(type=:af, N=3, pf=1);
ζgrid = init_grid(type=:af, N=1, pf=-1);

grids = init_grids(rgrid, θgrid, ζgrid);

solver = init_solver(full_spectrum=true, prob=prob);

evals, ϕft, ϕ = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=surfs);


#continuum_plot(evals)

tae_ind = find_ind(evals, 0.284)
tae_freq = evals.ω[tae_ind]

#potential_plot(ϕft, grids, tae_ind)
@test tae_ind == 28
@test tae_freq ≈ 0.284 atol=0.01

