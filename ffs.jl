
#ffs with an island seems to have no island modes...



Nr=30;
Nθ=10;


#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=1000.0);
isl = init_island(m0=2, n0=-1, w=0.1) #4e-4 is actually still pretty small.
prob = init_problem(q=island_mode_21, geo=geo, isl=isl)#, met=slab_metric!); #, met=no_delta_metric!)



rgrid = rfem_grid(N=Nr, sep1=0.42, sep2 =0.58, frac=0.4)
θgrid = afem_grid(N=Nθ, pf=2)
#ζgrid = afem_grid(N=Nζ, pf=-1)
ζgrid = asm_grid(start=-4, N=4)
grids = init_grids(rgrid, θgrid, ζgrid);


evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=false);

continuum_plot(evals)

ind = find_ind(evals, 0.1478122)

ind = 13

potential_plot(ϕft, grids, ind)


contour_plot(ϕ, grids, ind=ind)