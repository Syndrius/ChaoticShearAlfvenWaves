
sgrid = init_grid(:s, 15, start=0.15, stop=0.9)
ϑgrid = init_grid(:ϑ, 4, pf=1)
ζgrid = init_grid(:ζ, 1, pf=-1)

grids = init_grids(sgrid, ϑgrid, ζgrid)

solver = init_solver(full_spectrum=true, prob=prob)

evals, _, _ = compute_spectrum(prob, grids, solver, surfs);


tae_ind = find_ind(evals, 0.25)
tae_freq = evals.ω[tae_ind]


@test tae_ind == 28
@test tae_freq ≈ 0.24188 atol=0.001
