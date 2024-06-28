
#this file is fkd now, has been used to compare FSS to FFS.



#First benchmark, replicates the Simple TAE from Fu and Van Dam 1989
#they find a TAE with normalised frequency at Ω=0.31
#due to the coupling between (m, n)=(1, -1) and (2, -1)
#Note that they have Δ'=0 (No Shafranov shift)
using MID
using Plots; plotlyjs() #having this here, and installed in the global environment tricks it into using plotlyjs for interactive plots. This is an awful solution. Plotlyjs doesn't work with Latex strings, so can't really be used for final copies
#also gives some fkn warning, I think becuase MID doesn't have PlotlyJS.



N = 30;

geo = GeoParamsT(R0=4.0)

prob = init_problem(q=fu_dam_q, geo=geo, met=no_delta_metric!); 
#feels a little clunky, but I guess we can always create a init_grids function that takes the individual args later.
rgrid = init_fem_grid(N=N)
θgrid = init_sm_grid(start=1, count=2)
#θgrid = init_fem_grid(N=8, pf=1)
ζgrid = init_sm_grid(start=-1, count=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)
#tae_freq = (0.395845)^2/10^2

#inputs_to_file(prob=prob, grids=grids, dir="data/")

ω_cont = continuum(prob=prob, grids=grids);


plot_continuum(ω = ω_cont, grids=grids)


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true);

tae_ind = find_ind(ω, 0.31)
display(ω[tae_ind])




ϕms = mode_structure(ϕ, grids);
reconstruct_continuum(ω, ϕms, grids)#, filename="data/fu_dam_spectrum.png")
plot_potential(ϕms, grids, tae_ind, 1)

#this looks a little weird, can happen when the continuum bends.
#just an artifact of our simple reconstruction method.



reconstruct_continuum(ω, ϕ, grids)#, filename="data/fu_dam_spectrum.png")

plot_potential(ϕ, grids, tae_ind, 1)#, filename="data/fu_dam_tae.png")


#with R0=4, as per original case, we get 0.30999, or 0.31, 
#with R0=10, we get 0.3227 vs 0.3237
#with R0=20, we get 0.3276, vs 0.3283
#think we can be pretty confident.

#for Berks equation, ie with Δ,
#with R0=4, we get 0.3022, vs 0.2984
#with R0=10, we get 0.3188, vs 0.3178
#with R0=20, we get 0.3256, vs 0.3252

#again, pretty confident, Δ seems to shift tae freq down a bit, which is reduced as ϵ gest larger (smaller?? Bigger R0...)