
#TODO!

using MID
using Plots; plotlyjs()

#actually including Δ does seem to make our overestimated damping predictions more overestimated.
#something funny is going on here.

#with Δ'=0,
#ignore this, did a proper convergence test,
#very promising result just hasn't convereged yet. Seems to converge to ~0.0063 so a bit larger than other cases, as before.
#may need to run a proper convergence test with δ=0.
#ω = 0.32578522869815174 - 0.0055639292010700595im with N=4000, δ=-4e-11
#ω = 0.3268323919472274 - 0.006429238949748201im with N=5000, δ=-4e-11
#ω = 0.32604286776360225 - 0.0009537158133668805im with N=5000, δ=-4e-12
#ω = 0.3302203378605953 - 0.006120385872870546im with N=2000, δ=-4e-9
#ω = 0.32597136448865865 - 0.005571227966510881im with N=4000, δ=-4e-11, with only two modes!

#with test_metric 
#0.3255335613283627 - 0.006115903141937243im
#-0.018787319860296026

#normal metric
#0.32821226460709313 - 0.007118477121437926im
#-0.021688638387598165

N = 4000; 
#the collect is a bit annoying, but ok because we will typically use a clustered grid.
#rgrid = collect(LinRange(0, 1, N));
#rgrid = clustered_grid(N, 0.75, 0.8, 0.2)

geo = GeoParamsT(R0=10.0)

#isl = IslandT(A=4e-5, m0=5, n0=4)

#pretty confident it is convergeing to 0.32780760640103157 - 0.006921689729791451im
#giving ratio as -0.021115097986236567, so consistently above what the literature is giving!


prob = init_problem(q=singular_bowden_q, geo=geo, δ=-4.0e-9, dens=bowden_singular_dens)#, met=test_metric!); #probbaly should use geo if it is part 
#even if it is not really used.
grids = init_grids(N=N, sep1=0.75, sep2=0.8, frac=0.2, mstart=1, mcount=2, nstart=-1, ncount=1, f_quad=4);
#extra little spike is indep of clustered grid by the looks of it.
#maybe a concern but not sure what to do.
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1, f_quad=4);
#tae_freq = (0.381 / geo.R0)^2


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=(0.33/geo.R0)^2, reconstruct=true);

display(ω[1])
#display(abs(ω[1]))
display(imag(ω[1])/real(ω[1]))

reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)

tae_ind = find_ind(ω, 0.33)
tae_ind = 1
plot_potential(grids=grids, ϕ=ϕ, ind=tae_ind, n=1)

display(ω[tae_ind])
display(imag(ω[tae_ind])/real(ω[tae_ind]))

#converging to double somehow???
display(imag(ω[tae_ind]^2)/2)

tae_freq = (ω[tae_ind] / geo.R0)^2