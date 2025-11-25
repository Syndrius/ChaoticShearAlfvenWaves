
geo = init_geometry(:κ, R0=1.0)

isl = init_island(:κ, m0=1, n0=-1, w=0.1, ψ0=0.5, qp=1.0)

#need to change this to allow a single isl to be input.
fields = init_fields(:κ, q=island_q, isl=isl)
prob = init_problem(geometry=geo, fields=fields)

κgrid = init_grid(:κ, 20, start=0.001, stop=0.9999)
ᾱgrid = init_grid(:ᾱ, 5, pf=1)
τgrid = init_grid(:sm, 1, start=0)

grids = init_grids(κgrid, ᾱgrid, τgrid)
solver = init_solver(prob=prob, full_spectrum=true)
evals, _, _ = compute_spectrum(prob, grids, solver);

ind = find_ind(evals, 0.15)

@test ind == 77
@test real(evals.ω[ind]) ≈ 0.1424 atol=0.001

