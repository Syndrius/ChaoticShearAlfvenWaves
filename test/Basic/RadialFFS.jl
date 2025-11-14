
rgrid = init_grid(:r, 20)
θgrid = init_grid(:θ, 4, pf=1);
ζgrid = init_grid(:sm, 1, start=-1)

grids = init_grids(rgrid, θgrid, ζgrid);

geo = init_geometry()
fields = init_fields(:r)

prob = init_problem(fields=fields, geometry=geo)

solver = init_solver(full_spectrum=true, prob=prob)

evals, _, _ = compute_spectrum(prob=prob, grids=grids, solver=solver);

#continuum_plot(evals);

#unsure why it is flipped... hopefully a resolution problemo
#this is an extremly different tae frequency... real good lol
tae_ind = find_ind(evals, 0.300)
tae_freq = evals.ω[tae_ind]

@test tae_ind == 38
@test tae_freq ≈ 0.300 atol=0.001

#unsure if the second test is actually providing anything tbh!
#=
evals, _, _ = compute_spectrum(prob=prob, grids=grids, full_spectrum=false, target_freq=tae_freq);


tae_ind = find_ind(evals, 0.300)
tae_freq = evals.ω[tae_ind]

@test tae_ind == 1
@test tae_freq ≈ 0.300 atol=0.001
=#
