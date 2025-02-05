
#example with multiple islands

using MID


using Plots; plotlyjs()


#start very small, matrix scales much more extremly
Nr = 100;

geo = GeoParamsT(R0=10.0)


isl = IslandT(m0=2.0, n0=-1.0, A=0e-4)
isl2 = IslandT(m0=3.0, n0=-2.0, A=0e-4)
#rgrid = collect(LinRange(0, 1, N));
prob = init_problem(q=chaos_q, geo=geo, isl=isl, isl2=isl2); 
rgrid = rfem_grid(N=Nr)
θgrid = asm_grid(start=2, N = 2)
ζgrid = asm_grid(start=-2, N = 2)
grids = init_grids(rgrid, θgrid, ζgrid);
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
#grids = init_grids(rgrid, θgrid, ζgrid)



evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);

continuum_plot(evals)

tae_ind = find_ind(evals, 0.06327)

potential_plot(ϕft, grids, tae_ind)

contour_plot(ϕ, grids, ind=tae_ind)