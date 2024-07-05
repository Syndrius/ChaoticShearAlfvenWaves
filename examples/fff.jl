
using MID
using Plots; plotlyjs()

#Integration has been massivly sped up, but this is still slow af, probbaly requires multiproc.
#~90% of the time is spent numerically integrating. Wonder if there is anything we can do???
rgrid = init_fem_grid(N=20);
θgrid = init_fem_grid(N=5, pf=2);
ζgrid = init_fem_grid(N=3, pf=-2);

grids = init_grids(rgrid, θgrid, ζgrid);

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo); 


#with @views. 22.907080 seconds (8.07 M allocations: 721.379 MiB, 1.36% gc time)
#fk load more allocations and gc without views.
#outrageous spead up shifting the ϕ[:, test, :, ...] to ϕ[:, testr, testθ, :, :]

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true, reconstruct=true); 


ϕms = mode_structure(ϕ, grids);
reconstruct_continuum_n(ω, ϕms, grids, -2)#, ymax=0.25)#, filename="label_test.png")


#almost impossible to tell from the continuum reconstruction,
#but the same tae is found with both FFS and FFF which is v promising.
#even with tiny grid they seem to have v similar frequency and mode structure!
tae_ind = find_ind(ω, 0.396)
tae_freq = ω[tae_ind]
plot_potential(ϕms, grids, tae_ind, 1)