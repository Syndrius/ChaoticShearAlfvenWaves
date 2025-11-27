

geo = init_geometry(:tor, R0=10.0)
fields = init_fields(:r, q=damping_q, dens=damping_dens)
flr = init_flr(δ=-4.0e-9)
prob = init_problem(geometry=geo, fields=fields, flr=flr)

rgrid = init_grid(:r, 1000, sep1=0.75, sep2=0.78, frac=0.5)
θgrid = init_grid(:sm, 2, start=1)
ζgrid = init_grid(:sm, 1, start=-1)

grids = init_grids(rgrid, θgrid, ζgrid)

solver = init_solver(prob=prob, target=0.323, nev=5)
#%%
evals, _, _ = compute_spectrum(prob, grids, solver);

@test real(evals.ω[1]) ≈ 0.3231 atol=0.001
@test imag(evals.ω[1]) ≈ -0.004778 atol = 0.001
@test imag(evals.ω[1]) / real(evals.ω[1]) ≈ -0.014789 atol=0.001


