
Nr = 30;

geo = init_geo(R0=4.0)

prob = init_problem(q=fu_dam_q, geo=geo)#, met=no_delta_metric!); 

rgrid = init_grid(type=:rf, N=Nr)
θgrid = init_grid(type=:as, N=2, start=1)
ζgrid = init_grid(type=:as, N=1, start=-1)

grids = init_grids(rgrid, θgrid, ζgrid)

solver = init_solver(full_spectrum=true, prob=prob)

#evals, _, _ = compute_spectrum(prob=prob, grids=grids, solver=solver);

#not ideal
dir = "./"

inputs_to_file(dir=dir, prob=prob, solver=solver, grids=grids)

MID.Solve.spectrum_from_file(dir)

evals = evals_from_file(dir=dir);



#continuum_plot(evals)


tae_ind = find_ind(evals, 0.289)
tae_freq = evals.ω[tae_ind]

#potential_plot(ϕft, grids, tae_ind)

@test tae_ind == 58
@test tae_freq ≈ 0.2889 atol=0.001

rm(dir * "grids.jld2")
rm(dir * "prob.jld2")
rm(dir * "solver.jld2")
rm(dir * "evals.jld2")
