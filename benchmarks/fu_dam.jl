
#First benchmark, replicates the Simple TAE from Fu and Van Dam 1989
#they find a TAE with normalised frequency at Ω=0.31
#due to the coupling between (m, n)=(1, -1) and (2, -1)
#Note that they have Δ'=0 (No Shafranov shift)
using MID
#using Plots#; plotlyjs() #having this here, and installed in the global environment tricks it into using plotlyjs for interactive plots. This is an awful solution. Plotlyjs doesn't work with Latex strings, so can't really be used for final copies
#also gives some fkn warning, I think becuase MID doesn't have PlotlyJS.



N = 200;

geo = GeoParamsT(R0=4.0)

prob = init_problem(q=fu_dam_q, geo=geo, met=no_delta_metric!); 
grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
#tae_freq = (0.395845)^2/10^2

#inputs_to_file(prob=prob, grids=grids, dir="data/")

ω_cont = continuum(prob=prob, grids=grids);


plot_continuum(ω = ω_cont, grids=grids)


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true);

#this looks a little weird, can happen when the continuum bends.
#just an artifact of our simple reconstruction method.
reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)#, filename="data/fu_dam_spectrum.png")

tae_ind = find_ind(ω, 0.31)
display(ω[tae_ind])
plot_potential(grids=grids, ϕ=ϕ, ind=tae_ind, n=1)#, filename="data/fu_dam_tae.png")


#with R0=4, as per original case, we get 0.30999, or 0.31, 
#with R0=10, we get 0.3227 vs 0.3237
#with R0=20, we get 0.3276, vs 0.3283
#think we can be pretty confident.

#for Berks equation, ie with Δ,
#with R0=4, we get 0.3022, vs 0.2984
#with R0=10, we get 0.3188, vs 0.3178
#with R0=20, we get 0.3256, vs 0.3252

#again, pretty confident, Δ seems to shift tae freq down a bit, which is reduced as ϵ gest larger (smaller?? Bigger R0...)