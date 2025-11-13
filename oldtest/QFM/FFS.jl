
Nr = 20;
Nθ = 4

sgrid = init_grid(type=:rf, N=Nr, start=0.4, stop=0.9)
ϑgrid = init_grid(type=:af, N=Nθ, pf=1)
φgrid = init_grid(type=:as, N=1, start=-1)
grids = init_grids(sgrid, ϑgrid, φgrid)

solver = init_solver(full_spectrum=true, prob=prob)

evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=surfs);

#continuum_plot(evals);

#unsure why it is flipped... hopefully a resolution problemo
#this is an extremly different tae frequency... real good lol
tae_ind = find_ind(evals, 0.2815)
tae_freq = evals.ω[tae_ind]

#potential_plot(ϕft, grids, tae_ind)

@test tae_ind == 38
@test tae_freq ≈ 0.2815 atol=0.01
