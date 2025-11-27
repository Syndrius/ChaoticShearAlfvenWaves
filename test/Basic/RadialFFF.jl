

rgrid = init_grid(:r, 15)
θgrid = init_grid(:θ, 3, pf=1);
ζgrid = init_grid(:ζ, 1, pf=-1);

grids = init_grids(rgrid, θgrid, ζgrid);

solver = init_solver(full_spectrum=true, prob=prob);

#with @views. 22.907080 seconds (8.07 M allocations: 721.379 MiB, 1.36% gc time)
#fk load more allocations and gc without views.
#outrageous spead up shifting the ϕ[:, test, :, ...] to ϕ[:, testr, testθ, :, :]

evals, _, _ = compute_spectrum(prob, grids, solver);


#continuum_plot(evals)


#almost impossible to tell from the continuum reconstruction,
#but the same tae is found with both FFS and FFF which is v promising.
#even with tiny grid they seem to have v similar frequency and mode structure!
tae_ind = find_ind(evals, 0.301)
tae_freq = evals.ω[tae_ind]

@test tae_ind == 28
@test tae_freq ≈ 0.301 atol=0.001
