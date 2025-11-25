
#testing with flux
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

#testing with geometric radius.
rgrid = init_grid(:r, 15)
θgrid = init_grid(:θ, 3, pf=1);
ζgrid = init_grid(:ζ, 1, pf=-1);

grids = init_grids(rgrid, θgrid, ζgrid);

geo = init_geometry()
fields = init_fields(:r)

prob = init_problem(fields=fields, geometry=geo)

solver = init_solver(full_spectrum=true, prob=prob);
solver = init_solver(prob=prob, targets=[0.301, 0.4], nev=5)

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

#@test tae_ind == 28
@test tae_freq ≈ 0.301 atol=0.001
