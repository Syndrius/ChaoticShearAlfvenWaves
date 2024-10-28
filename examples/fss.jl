
using MID
using Plots; plotlyjs()


#start very small, matrix scales much more extremly
Nr = 100;

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo)#, met=no_delta_metric!); 

rgrid = rfem_grid(N=Nr)
θgrid = asm_grid(start=2, N=2)
ζgrid = asm_grid(start=-2, N=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)



evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);

#scatter(cr.r, real.(cr.ω), ylimits=(-0.05, 1.05))
continuum_plot(evals)

ind = find_ind(evals, 0.3762)
ind = 348
plot_potential(ϕft, grids, ind)


contour_plot(ϕ, grids, ind=ind)




reconstruct_continuum(ω, ϕ, grids)#, filename="data/fu_dam_spectrum.png")

#unsure why it is flipped... hopefully a resolution problemo
tae_ind = find_ind(ω, 0.39)
tae_freq = ω[tae_ind]
plot_potential(ϕ, grids, tae_ind, 1)


#now increase the resolution,
Nr = 500

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo)#, met=no_delta_metric!); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=2, count=2)
ζgrid = init_sm_grid(start=-2, count=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)

#tae_freq=0.30
ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=0.376);



reconstruct_continuum(ω, ϕ, grids)#, filename="data/fu_dam_spectrum.png")

#unsure why it is flipped... hopefully a resolution problemo
tae_ind = find_ind(ω, 0.38)
tae_freq = ω[tae_ind]
#now tae modes are same sign, is that correct??? would match the fem1d case.

#the mode structure for the tae look identical, but the surface is wildly different
#mode structure for non-tae modes is much more localized, is that good?
#why are these so different?
#non-tae modes seem to look better, but idk wth.
plot_potential(ϕ, grids, tae_ind, 1)

ϕsur = construct_surface(ϕ, length(ω), grids, 0.0);


#this one looks more like a combo of m=2 and m=3 but fk me 
#we must be doing something wrong with the fourier transforms!
plot_phi_surface(ϕsur, grids, tae_ind)