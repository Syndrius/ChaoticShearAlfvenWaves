
Nr = 20;
Nθ = 4


rgrid = init_grid(type=:rf, N=Nr, start=0.4)
θgrid = init_grid(type=:af, N=Nθ, pf=1)
ζgrid = init_grid(type=:as, N=1, start=-1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)

solver = init_solver(full_spectrum=true, prob=prob)

evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=surfs);

#continuum_plot(evals);

#unsure why it is flipped... hopefully a resolution problemo
#this is an extremly different tae frequency... real good lol
tae_ind = find_ind(evals, 0.284)
tae_freq = evals.ω[tae_ind]

#potential_plot(ϕft, grids, tae_ind)

@test tae_ind == 38
@test tae_freq ≈ 0.284 atol=0.01
