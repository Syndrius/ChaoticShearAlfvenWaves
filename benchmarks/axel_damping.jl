

#file for comparing the damping computed with MID, vs with the direct solution of Berks simplified form.
#TODO!

using MID
using Plots; plotlyjs()

#with og W
#0.3911009635290564 - 0.0019341830644314319im
#-0.001512921720281432

#full metric, N=1000, δ=7
#0.39266476207502504 - 0.0021696044043428894im

#-0.0017038543944564541

#this is giving ~1.8 times what we expect.


#full met, N=3000, δ=9, R0=10
#0.39265080894305077 - 0.0021263116078717im
#-0.001669795945791644
#R0=20
#0.3960742171412108 - 0.0007261131909515339im
#-0.0005751894273240706

#test_met, N=3000, δ=9, R0=10
#0.36420118781247635 - 0.0014092290709537005im
#-0.0010264858030824204
#R0=20

#increasing to R0=20 and comparing the the Axel contour method in phi_two_mode
#our code with diagonal metrix is v v close. %diff of ~1.7%

N = 1000; 
#the collect is a bit annoying, but ok because we will typically use a clustered grid.
#rgrid = collect(LinRange(0, 1, N));

#getting some weird spikes, may need to double check this is working as expected??
#rgrid = clustered_grid(N, 0.91, 0.98, 0.25)

geo = GeoParamsT(R0=10.0)

#isl = IslandT(A=4e-5, m0=5, n0=4)

#pretty confident it is convergeing to 0.32780760640103157 - 0.006921689729791451im
#giving ratio as -0.021115097986236567, so consistently above what the literature is giving!
#test metric is pretty cooked, real tae freq is pretty close, touch higher with og metric
#damping is just completly cooked though, may need to check out damping 
#implementation.
prob = init_problem(q=Axel_q, geo=geo, δ=-4.0e-7, dens=axel_dens, met=diagonal_toroidal_metric!); #probbaly should use geo if it is part of prob,
#prob = init_problem(q=singular_bowden_q, geo=geo, δ=-4e-9, dens=bowden_singular_dens); #probbaly should use geo if it is part 
#even if it is not really used.
#grids = init_grids(N=N, sep1=0.91, sep2=0.98, frac=0.25, mstart=2, mcount=2, nstart=-2, ncount=1);
rgrid = init_fem_grid(N=N)
θgrid = init_sm_grid(start=2, count=2)
ζgrid = init_sm_grid(start=-2, count=1)
grids = init_grids(rgrid, θgrid, ζgrid)
#grids = init_grids(N=N, mstart=2, mcount=2, nstart=-2, ncount=1);
#tae_freq = (0.381 / geo.R0)^2


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=(0.390/geo.R0)^2, reconstruct=true);
tae_ind = 1
display(ω[tae_ind])

display(imag(ω[tae_ind]^2))



#scale of eigenfunctions is v different here.
#also singularity is much much more pronounced...

reconstruct_continuum(ω = ω.^2, ϕ = ϕ, grids = grids, ymax=0.25, ymin=-0.02)

tae_ind = find_ind(ω.^2, 0.157)

tae_ind = 1
plot_potential(grids=grids, ϕ=ϕ, ind=tae_ind, n=1)

display(ω[tae_ind])
display(imag(ω[tae_ind])/real(ω[tae_ind]))

#converging to double somehow???
display(imag(ω[tae_ind]^2)/2)

tae_freq = (ω[tae_ind] / geo.R0)^2



ω, ϕ = analytical_construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=(0.395/geo.R0)^2, reconstruct=true);

tae_ind = 1
plot_potential(grids=grids, ϕ=ϕ, ind=tae_ind, n=1)

display(ω[tae_ind])

display(imag(ω[tae_ind]^2))