Nr = 15;
Nθ = 3
Nζ = 1

sgrid = init_grid(type=:rf, N=Nr, start=0.4, stop=0.9)
ϑgrid = init_grid(type=:af, N=Nθ, pf=1)
φgrid = init_grid(type=:af, N=Nζ, pf=-1)
grids = init_grids(sgrid, ϑgrid, φgrid)

solver = init_solver(full_spectrum=true, prob=prob);

evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=surfs);


#continuum_plot(evals)

tae_ind = find_ind(evals, 0.2818)
tae_freq = evals.ω[tae_ind]

#potential_plot(ϕft, grids, tae_ind)
@test tae_ind == 28
@test tae_freq ≈ 0.2818 atol=0.01

