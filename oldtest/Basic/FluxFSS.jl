

#start very small, matrix scales much more extremly
Nψ = 30;

geo = init_geo(R0=4.0)

prob = init_problem(type=:flux, q=MID.Equilibrium.flux_fu_dam_q, geo=geo)#, met=no_delta_metric!); 

ψgrid = init_grid(type=:rf, N=Nψ, stop=0.5)
θgrid = init_grid(type=:as, N=2, start=1)
ζgrid = init_grid(type=:as, N=1, start=-1)

grids = init_grids(ψgrid, θgrid, ζgrid)

solver = init_solver(full_spectrum=true, prob=prob)

evals, _, _ = compute_spectrum(prob=prob, grids=grids, solver=solver);

#continuum_plot(evals)


tae_ind = find_ind(evals, 0.291)
tae_freq = evals.ω[tae_ind]

#potential_plot(ϕft, grids, tae_ind)

@test tae_ind == 58
@test tae_freq ≈ 0.291 atol=0.001


#=
evals, _, _ = compute_spectrum(prob=prob, grids=grids, full_spectrum=false, target_freq = tae_freq);


tae_ind = find_ind(evals, 0.2879)
tae_freq = evals.ω[tae_ind]

@test tae_ind == 1
@test tae_freq ≈ 0.2879 atol=0.001
=#
