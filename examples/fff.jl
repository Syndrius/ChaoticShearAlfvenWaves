
using MID
using Plots; plotlyjs()

#Integration has been massivly sped up, but this is still slow af, probbaly requires multiproc.
#~90% of the time is spent numerically integrating. Wonder if there is anything we can do???
Nr = 30
Nθ = 6
Nζ = 1
rgrid = init_fem_grid(N=Nr);
θgrid = init_fem_grid(N=Nθ, pf=2);
ζgrid = init_fem_grid(N=Nζ, pf=-2);

grids = init_grids(rgrid, θgrid, ζgrid);

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo); 


#with @views. 22.907080 seconds (8.07 M allocations: 721.379 MiB, 1.36% gc time)
#fk load more allocations and gc without views.
#outrageous spead up shifting the ϕ[:, test, :, ...] to ϕ[:, testr, testθ, :, :]

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true); 

plot_continuum(evals)


ind = find_ind(evals, 0.383)
#ind = 348
plot_potential(ϕft, grids, ind)


contour_plot(ϕ, grids, ind=ind)