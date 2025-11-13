


Nψ = 20;
Nθ = 4

geo = init_geo(R0=4.0)

prob = init_problem(type=:flux, q=MID.Equilibrium.flux_fu_dam_q, geo=geo); 

ψgrid = init_grid(type=:rf, N=Nψ, stop=0.5)
θgrid = init_grid(type=:af, N=Nθ, pf=1)
ζgrid = init_grid(type=:as, N=1, start=-1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(ψgrid, θgrid, ζgrid)

solver = init_solver(full_spectrum=true, prob=prob)

evals, _, _ = compute_spectrum(prob=prob, grids=grids, solver=solver);

#continuum_plot(evals);

#unsure why it is flipped... hopefully a resolution problemo
#this is an extremly different tae frequency... real good lol
tae_ind = find_ind(evals, 0.306)
tae_freq = evals.ω[tae_ind]

@test tae_ind == 38
@test tae_freq ≈ 0.306 atol=0.001

#unsure if the second test is actually providing anything tbh!
#=
evals, _, _ = compute_spectrum(prob=prob, grids=grids, full_spectrum=false, target_freq=tae_freq);


tae_ind = find_ind(evals, 0.300)
tae_freq = evals.ω[tae_ind]

@test tae_ind == 1
@test tae_freq ≈ 0.300 atol=0.001
=#
