
#see how fss does with the singular case, not going for damping yet lol, that is to hard..

using MID
using Plots; plotlyjs()


#start very small, matrix scales much more extremly
Nr = 30;
Nθ = 7

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=bowden_singular_q, geo=geo, dens=bowden_singular_dens); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=1)
ζgrid = init_sm_grid(start=-1, count=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)



ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true);


ϕms = mode_structure(ϕ, grids);
reconstruct_continuum(ω, ϕms, grids)#, filename="data/fu_dam_spectrum.png")

#unsure why it is flipped... hopefully a resolution problemo
#this is an extremly different tae frequency... real good lol
tae_ind = find_ind(ω, 0.32)
tae_freq = ω[tae_ind]
plot_potential(ϕms, grids, tae_ind, 1)


#now increase the resolution,
#note with Nr=300, Nθ=30, this is almost a perfect prediction!
#including damping lol! -> v surprising since this is quite a low res!
#500x10 also gives a pretty good estimate. Tis still a bit off, we will need to move to parallel 
#to properly check this.
Nr = 500
Nθ = 10

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=bowden_singular_q, geo=geo, dens=bowden_singular_dens, met=diagonal_toroidal_metric!, δ=-4.0e-7);  

rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=1)
ζgrid = init_sm_grid(start=-1, count=1)
grids = init_grids(rgrid, θgrid, ζgrid)


tae_freq = 0.325

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq, reconstruct=true);


ϕms = mode_structure(ϕ, grids);
reconstruct_continuum(ω, ϕms, grids)#, filename="data/fu_dam_spectrum.png")

#unsure why it is flipped... hopefully a resolution problemo
tae_ind = find_ind(ω, 0.325)
tae_ind = 1
tae_freq = ω[tae_ind]

display(imag(ω[tae_ind])/real(ω[tae_ind]))
#now tae modes are same sign, is that correct??? would match the fem1d case.

#the mode structure for the tae look identical, but the surface is wildly different
#mode structure for non-tae modes is much more localized, is that good?
#surely the problem has to be with the phase factor stuff?
#doens't explain why the frequency is so cooked though
#maybe we can replicate this case with cka in 3d?? should be small enough right?
#why are these so different?
#trying with pf=0 is hard to understand... not sure what is going on tbh.
#pf = 1 seems somewhat in between the pf=2 case and the fss case. but much closer to pf=2.
tae_ind = 1
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