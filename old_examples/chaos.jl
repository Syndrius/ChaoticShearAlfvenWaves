
#example with multiple islands

using MID


using Plots; plotlyjs()


#start very small, matrix scales much more extremly
#Nr = 100;

geo = GeoParamsT(R0=10.0)


isl = IslandT(m0=2.0, n0=-1.0, A=0.0e-5)
isl2 = IslandT(m0=3.0, n0=-2.0, A=0.0e-5)
#rgrid = collect(LinRange(0, 1, N));
prob = init_problem(q=chaos_q, geo=geo, isl=isl, isl2=isl2); 
#rgrid = rfem_grid(N=Nr)
#θgrid = asm_grid(start=2, N = 2)
#ζgrid = asm_grid(start=-2, N = 2)
#grids = init_grids(rgrid, θgrid, ζgrid);

Nr = 30
Nθ = 6
Nζ = 2
rgrid = rfem_grid(N=Nr, gp=4, sep1=0.4, sep2=0.6, frac=0.0);
θgrid = afem_grid(N=Nθ, pf=2, gp=4);
ζgrid = afem_grid(N=Nζ, pf=-2, gp=4);
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
#grids = init_grids(rgrid, θgrid, ζgrid)



evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=false, target_freq=0.3);

continuum_plot(evals)

tae_ind = find_ind(evals, 0.06327)

potential_plot(ϕft, grids, tae_ind)

contour_plot(ϕ, grids, ind=tae_ind)

#haveing a look at the og continuum

Nr = 100;
geo = GeoParamsT(R0=10000.0)

prob = init_problem(q=chaos_q, geo=geo)#, met=no_delta_metric!); 

#rgrid = init_fem_grid(N=Nr)
rgrid = MID.ContGridDataT(N=100)
θgrid = asm_grid(start=-5, N=11)
ζgrid = asm_grid(start=-2, N=5)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)


ω_cont = continuum(prob, grids);

continuum_plot(ω_cont, grids, ymax=1)