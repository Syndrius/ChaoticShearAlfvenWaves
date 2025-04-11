
rgrid = MID.Structures.rfem_grid(N=30, start=0.5)#, stop=0.7)
θgrid = MID.Structures.asm_grid(start=1, N=2)#, f_quad=1)
ζgrid = MID.Structures.asm_grid(start=-1, N=1)#, f_quad=1)#, incr=2)

grids = init_grids(rgrid, θgrid, ζgrid)

#%%

#solver = init_solver(nev=100, target = 0.3, prob=prob)
solver = init_solver(full_spectrum=true, prob=prob)

#%%

evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=surfs);

#continuum_plot(evals)

tae_ind = find_ind(evals, 0.277)
tae_freq = evals.ω[tae_ind]

#potential_plot(ϕft, grids, tae_ind)

@test tae_ind == 58
@test tae_freq ≈ 0.277 atol=0.01

