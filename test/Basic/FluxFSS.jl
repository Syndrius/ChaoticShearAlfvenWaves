

ψgrid = init_grid(:ψ, 30)
θgrid = init_grid(:sm, 2, start=1)
ζgrid = init_grid(:sm, 1, start=-1)


grids = init_grids(ψgrid, θgrid, ζgrid)

solver = init_solver(full_spectrum=true, prob=prob)

evals, _, _ = compute_spectrum(prob, grids, solver);


tae_ind = find_ind(evals, 0.2717)
tae_freq = evals.ω[tae_ind]


@test tae_ind == 58
@test tae_freq ≈ 0.2717 atol=0.001

