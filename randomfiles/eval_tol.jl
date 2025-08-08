
using Arpack
using MID
#%%

#first define the problem
geo = init_geo(R0=4.0)

#rgrid = collect(LinRange(0, 1, N));
prob = init_problem(q=fu_dam_q, geo=geo); 
#%%
#then create the grids
Nr = 30;
rgrid = init_grid(type=:rf, N=Nr)
θgrid = init_grid(type=:as, N = 2, start = 1)
ζgrid = init_grid(type=:as, N = 1, start = -1)
grids = init_grids(rgrid, θgrid, ζgrid);

#%%
#then define the solver
#solver = init_solver(full_spectrum=true, prob=prob)
#solver = init_solver(target=0.33, prob=prob)
#solver = init_solver(prob=prob, nshifts=4)
solver = init_solver(prob=prob, targets=[0.2, 0.21], nev=20)
#%%


W, I = MID.Construct.construct(prob, grids)

ev, ef, nconv, niter, nmult, resid = eigs(W, I, nev=solver.nev, sigma=solver.targets[2], tol=1e-10)

display(ev)
display(resid[1:10])
display(nconv)





evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);
#%%
#scatter(cr.r, real.(cr.ω), ylimits=(-0.05, 1.05))
continuum_plot(evals)
display(evals.ω)
