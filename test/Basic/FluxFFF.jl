

geo = init_geometry()
fields = init_fields()


prob = init_problem(fields=fields, geometry=geo)

ψgrid = init_grid(:ψ, 15)
θgrid = init_grid(:θ, 4, pf=1)
ζgrid = init_grid(:ζ, 1, pf=-1)

grids = init_grids(ψgrid, θgrid, ζgrid)

solver = init_solver(full_spectrum=true, prob=prob)

evals, _, _ = compute_spectrum(prob, grids, solver);


tae_ind = find_ind(evals, 0.2709)
tae_freq = evals.ω[tae_ind]


@test tae_ind == 28
@test tae_freq ≈ 0.2679 atol=0.001


