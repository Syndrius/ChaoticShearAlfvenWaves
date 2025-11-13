#test running the code from file
#this is a closer representation to how the parallel version works
using MID
using MIDViz
#%%
#first define the problem
geo = init_geo(R0=4.0)

#rgrid = collect(LinRange(0, 1, N));
prob = init_problem(q=fu_dam_q, geo=geo); 
#%%
#then create the grids
Nr = 100;
rgrid = init_grid(type=:rf, N=Nr)
θgrid = init_grid(type=:as, N = 2, start = 1)
ζgrid = init_grid(type=:as, N = 1, start = -1)
grids = init_grids(rgrid, θgrid, ζgrid);

#%%
#then define the solver
#solver = init_solver(full_spectrum=true, prob=prob)
#solver = init_solver(target=0.33, prob=prob)
#solver = init_solver(prob=prob, nshifts=4)
solver = init_solver(prob=prob, targets=[0.2, 0.3, 0.4, 0.5])

#%%
dir = "/Users/matt/phd/MID/data/"

inputs_to_file(dir=dir, solver=solver, prob=prob, grids=grids)

MID.Solve.spectrum_from_file(dir)

evals = evals_from_file(dir=dir);

continuum_plot(evals)
