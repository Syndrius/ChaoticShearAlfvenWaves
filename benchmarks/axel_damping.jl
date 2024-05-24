

#file for comparing the damping computed with MID, vs with the direct solution of Berks simplified form.
#TODO!

using MID
using Plots; plotlyjs()

#this is giving ~1.8 times what we expect.

N = 3000; 
#the collect is a bit annoying, but ok because we will typically use a clustered grid.
#rgrid = collect(LinRange(0, 1, N));

#getting some weird spikes, may need to double check this is working as expected??
rgrid = clustered_grid(N, 0.91, 0.98, 0.25)

geo = GeoParamsT(R0=10.0)

#isl = IslandT(A=4e-5, m0=5, n0=4)

#pretty confident it is convergeing to 0.32780760640103157 - 0.006921689729791451im
#giving ratio as -0.021115097986236567, so consistently above what the literature is giving!

prob = init_problem(q=Axel_q, geo=geo, δ=-4.0e-9, dens=axel_dens); #probbaly should use geo if it is part of prob,
#prob = init_problem(q=singular_bowden_q, geo=geo, δ=-4e-9, dens=bowden_singular_dens); #probbaly should use geo if it is part 
#even if it is not really used.
grids = init_grids(rgrid=rgrid, mstart=1, mcount=4, nstart=-2, ncount=1);
#tae_freq = (0.381 / geo.R0)^2


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=(0.395/geo.R0)^2, reconstruct=true);

display(ω[tae_ind])

display(imag(ω[tae_ind]^2))

reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)

tae_ind = find_ind(ω, 0.398)
tae_ind = 1
plot_potential(r=rgrid, ϕ=ϕ, ind=tae_ind, pmd=grids.pmd, n=1)

display(ω[tae_ind])
display(imag(ω[tae_ind])/real(ω[tae_ind]))

#converging to double somehow???
display(imag(ω[tae_ind]^2)/2)

tae_freq = (ω[tae_ind] / geo.R0)^2