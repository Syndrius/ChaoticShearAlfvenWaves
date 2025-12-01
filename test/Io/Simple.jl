
geo = init_geometry()
fields = init_fields()

prob = init_problem(fields=fields, geometry=geo)

ψgrid = init_grid(:ψ, 30)
θgrid = init_grid(:sm, 2, start=1)
ζgrid = init_grid(:sm, 1, start=-1)


grids = init_grids(ψgrid, θgrid, ζgrid)

solver = init_solver(full_spectrum=true, prob=prob)

dir = abspath(joinpath(pathof(ChaoticShearAlfvenWaves), "../../test/data/"))

inputs_to_file(dir=dir, grids=grids, prob=prob, solver=solver)

compute_spectrum(dir);

evals = evals_from_file(dir)

tae_ind = find_ind(evals, 0.2717)
tae_freq = evals.ω[tae_ind]


@test tae_ind == 58
@test tae_freq ≈ 0.2717 atol=0.001

rm(joinpath(dir, "evals.jld2"))
rm(joinpath(dir, "grids.jld2"))
rm(joinpath(dir, "solver.jld2"))
rm(joinpath(dir, "prob.jld2"))
rm(joinpath(dir, "efuncs/"), recursive=true)
rm(joinpath(dir, "efuncs_ft/"), recursive=true)

