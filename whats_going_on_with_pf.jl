
using MID


using Plots; plotlyjs()

#this file has showed that pf is cooked
#and has allowed plotting to be fixed
#we have also fixed periodicity
#so our efuncs are all looking the same now.

#start very small, matrix scales much more extremly
Nr = 100;

geo = GeoParamsT(R0=4.0)

prob = init_problem(q=fu_dam_q, geo=geo)#, met=no_delta_metric!); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=0, count=6)
ζgrid = init_sm_grid(start=-1, count=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids_fss = init_grids(rgrid, θgrid, ζgrid)



evals_fss, ϕ_fss, ϕft_fss = compute_spectrum(prob=prob, grids=grids_fss, full_spectrum=true);

plot_continuum(evals_fss)

ind = find_ind(evals_fss, 0.82537)

plot_potential(ϕft_fss, grids_fss, ind)
contour_plot(ϕ_fss, grids_fss, ind=ind)
surface_plot(ϕ_fss, grids_fss, ind=ind)

tae_ind_fss = find_ind(evals_fss, 0.29)


plot_potential(ϕft_fss, grids_fss, tae_ind_fss)

#well I'll be a monkey's uncle.
#using ifft to recreate ϕ actually gives the same fkn thing.
#that is heaps cooked bro.
#the two modes are pi/4 out of phase, which is kind of wild.
#still need to figure out wot the fk is going on with pf though.
#also this version does not handle periodicity properly.
#perhaps we should incldue that in the post-processing step.
#well that is good.
contour_plot(ϕ_fss, grids_fss, ind=tae_ind_fss)


#now compare with ffs case

Nr = 200;
Nθ = 30
geo = GeoParamsT(R0=4.0)
prob = init_problem(q=fu_dam_q, geo=geo); 

rgrid = init_fem_grid(N=Nr)
#so pf is cooked af yay.
θgrid = init_fem_grid(N=Nθ)#, pf=1)
ζgrid = init_sm_grid(start=-1, count=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids_ffs = init_grids(rgrid, θgrid, ζgrid)

#ok so the pf is actually cooking the mode structure....
#that is pretty fkn annoying. I guess we will have to ignore for now.
#not sure how to fix...

evals_ffs, ϕ_ffs, ϕft_ffs = compute_spectrum(prob=prob, grids=grids_ffs, full_spectrum=false, target_freq=0.31);

plot_continuum(evals_ffs)


ind = find_ind(evals_ffs, 0.218343)

plot_potential(ϕft_ffs, grids_ffs, ind)
contour_plot(ϕ_ffs, grids_ffs, ind=ind)
surface_plot(ϕ_ffs, grids_ffs, ind=ind)



tae_ind_ffs = find_ind(evals_ffs, 0.31)

plot_potential(ϕft_ffs, grids_ffs, tae_ind_ffs)


contour_plot(ϕ_ffs, grids_ffs, ind=tae_ind_ffs)

surface_plot(ϕ_ffs, grids_ffs, ind=tae_ind_ffs)




Nr = 50;
Nθ = 10
Nζ = 2
geo = GeoParamsT(R0=4.0)
prob = init_problem(q=fu_dam_q, geo=geo); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ)
ζgrid = init_fem_grid(N=Nζ, pf=-1) #need this despit it not working
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids_fff = init_grids(rgrid, θgrid, ζgrid)

evals_fff, ϕ_fff, ϕft_fff = compute_spectrum(prob=prob, grids=grids_fff, full_spectrum=false, target_freq = 0.31);

plot_continuum(evals_fff)


ind = find_ind(evals_fff, 0.3140)

plot_potential(ϕft_fff, grids_fff, ind)
contour_plot(ϕ_fff, grids_fff, ind=ind)
surface_plot(ϕ_fff, grids_fff, ind=ind)

tae_ind_fff = find_ind(evals_fff, 0.3011)
#fk me this is perf now!
#this should be our basic example now
plot_potential(ϕft_fff, grids_fff, tae_ind_fff)

#so I guess we don't actually know what a tae looks like or should look like?
#why the fk does this look like the old one...
#higher res fixed this.
contour_plot(ϕ_fff, grids_fff, ind=tae_ind_fff)
surface_plot(ϕ_fff, grids_fff, ind=tae_ind_fff)