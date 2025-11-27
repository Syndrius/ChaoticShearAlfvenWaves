
rgrid = init_grid(:r, 30)
θgrid = init_grid(:sm, 2, start=1)
ζgrid = init_grid(:sm, 1, start=-1)

grids = init_grids(rgrid, θgrid, ζgrid);

solver = init_solver(full_spectrum=true, prob=prob)

evals, _, _ = compute_spectrum(prob, grids, solver);


tae_ind = find_ind(evals, 0.289)
tae_freq = evals.ω[tae_ind]

#potential_plot(ϕft, grids, tae_ind)

@test tae_ind == 58
@test tae_freq ≈ 0.2889 atol=0.001


#=
evals, _, _ = compute_spectrum(prob=prob, grids=grids, full_spectrum=false, target_freq = tae_freq);


tae_ind = find_ind(evals, 0.2879)
tae_freq = evals.ω[tae_ind]

@test tae_ind == 1
@test tae_freq ≈ 0.2879 atol=0.001
=#
