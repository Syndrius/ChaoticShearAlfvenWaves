
#test with island anonymous functions.
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

dir = abspath(joinpath(pathof(ChaoticShearAlfvenWaves), "../../test/data/"))

inputs_to_file(dir=dir, grids=grids, prob=prob, solver=solver)

compute_spectrum(dir);

evals = evals_from_file(dir)

ind = find_ind(evals, 0.015)

@test ind == 77
@test real(evals.ω[ind]) ≈ 0.01424 atol=0.001

rm(joinpath(dir, "evals.jld2"))
rm(joinpath(dir, "grids.jld2"))
rm(joinpath(dir, "solver.jld2"))
rm(joinpath(dir, "prob.jld2"))
rm(joinpath(dir, "efuncs/"), recursive=true)
rm(joinpath(dir, "efuncs_ft/"), recursive=true)

