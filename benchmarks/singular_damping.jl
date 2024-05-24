
#TODO!

using MID
using Plots; plotlyjs()

#should split this into 3, so we can stop changing shit all the time.

#think the comparison method probably shows us that there is a problem, we are pretty consistently getting the wrong tae freq and overestimating the damping by ~2.
#may be time to tackle the weak form...

#if this is the converged value 0.3258893681504531 - 0.005535297205358447im
#that is probably close enough..., this required N=4000, m=0-3, and δ=-4e-11.
#this was done with the new weak form W. 
#given Axel case is still not working, it is unlikley that this is the true converged result..

N = 4000; 
#the collect is a bit annoying, but ok because we will typically use a clustered grid.
#rgrid = collect(LinRange(0, 1, N));
rgrid = clustered_grid(N, 0.75, 0.8, 0.2)

geo = GeoParamsT(R0=10.0)

#isl = IslandT(A=4e-5, m0=5, n0=4)

#pretty confident it is convergeing to 0.32780760640103157 - 0.006921689729791451im
#giving ratio as -0.021115097986236567, so consistently above what the literature is giving!


prob = init_problem(q=singular_bowden_q, geo=geo, δ=-4e-11, dens=bowden_singular_dens); #probbaly should use geo if it is part 
#even if it is not really used.
grids = init_grids(rgrid=rgrid, mstart=0, mcount=4, nstart=-1, ncount=1);
#tae_freq = (0.381 / geo.R0)^2


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=(0.3259/geo.R0)^2, reconstruct=true);

display(ω[tae_ind])
display(abs(ω[1]))
display(imag(ω[1])/real(ω[1]))

reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)

tae_ind = find_ind(ω, 0.295)
tae_ind = 1
plot_potential(r=rgrid, ϕ=ϕ, ind=tae_ind, pmd=grids.pmd, n=1)

display(ω[tae_ind])
display(imag(ω[tae_ind])/real(ω[tae_ind]))

#converging to double somehow???
display(imag(ω[tae_ind]^2)/2)

tae_freq = (ω[tae_ind] / geo.R0)^2