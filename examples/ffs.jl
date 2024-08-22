
using MID
using Plots; plotlyjs()

#with individual parts.
#different ways of doing it are basically identical.
#10.992975 seconds (15.38 M allocations: 1.195 GiB, 4.40% gc time)

#start very small, matrix scales much more extremly
Nr = 40;
Nθ = 6

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=2)
ζgrid = init_sm_grid(start=-2, count=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)



evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);

#scatter(cr.r, real.(cr.ω), ylimits=(-0.05, 1.05))
plot_continuum(evals)

#ϕms = mode_structure(ϕ, grids);
#reconstruct_continuum(ω, ϕms, grids)#, ymax=3)#, filename="data/fu_dam_spectrum.png")

#unsure why it is flipped... hopefully a resolution problemo
#this is an extremly different tae frequency... real good lol
tae_ind = find_ind(evals, 0.396)
tae_freq = evals.ω[tae_ind]

#ind = 7
plot_potential(ϕft, grids, tae_ind)

contour_plot(ϕ, grids, tae_ind)
#so not symmetric or the correct mode number????? wot.
#plot_phi_surface(ϕ, grids, ind)

#phase factor seems to be fine, need to check with island modes though!



#now increase the resolution,
Nr = 50
Nθ = 5

geo = GeoParamsT(R0=10.0)
prob = init_problem(q=Axel_q, geo=geo)#, met=no_delta_metric!); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=2)
ζgrid = init_sm_grid(start=-2, count=1)
grids = init_grids(rgrid, θgrid, ζgrid)


tae_freq = 0.39

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq, reconstruct=true);


ϕms = mode_structure(ϕ, grids);
reconstruct_continuum(ω, ϕms, grids)#, filename="data/fu_dam_spectrum.png")

#unsure why it is flipped... hopefully a resolution problemo
tae_ind = find_ind(ω, 0.396)
tae_freq = ω[tae_ind]
#now tae modes are same sign, is that correct??? would match the fem1d case.

#the mode structure for the tae look identical, but the surface is wildly different
#mode structure for non-tae modes is much more localized, is that good?
#surely the problem has to be with the phase factor stuff?
#doens't explain why the frequency is so cooked though
#maybe we can replicate this case with cka in 3d?? should be small enough right?
#why are these so different?
#trying with pf=0 is hard to understand... not sure what is going on tbh.
#pf = 1 seems somewhat in between the pf=2 case and the fss case. but much closer to pf=2.

#non-tae modes seem to look better, but idk wth.
plot_potential(ϕms, grids, tae_ind, 1)

#some pf scaling seems to have made the surface even less periodic which is kinda wild.
#seems to be a bit better now, which is good!
#freq is still completly cooked though!

#m labels are wrong, something to do with either pf or becuase it should start at 0? maybe both maybe neither.

#in fu-dam case one of the modes is much less prevalent that I think it should be!

#tae seems to approximatly match fss method by eye test, does look `less` periodic though!
#still need to understand the huge frequency difference.

#think this needs to factor the phase factor stuff, unclear how though!
plot_phi_surface(ϕ, grids, tae_ind, 1)