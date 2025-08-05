
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

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);
#%%
#scatter(cr.r, real.(cr.ω), ylimits=(-0.05, 1.05))
continuum_plot(evals)


ind = find_ind(evals, 0.289)
#ind = 348
potential_plot(ϕft, grids, ind)

#%%
C = [[1, 2]  [3, 4]  [5, 6]]
D = [[0.1, 0.2] [0.3, 0.4]]

T = C' * D
res = T * C
#%%
W = zeros(3, 3)

for i in 1:3, j in 1:3, k in 1:2, l in 1:2
    W[i, j] += C[k, i] * D[k, l] * C[l, j]
end
display(W)
