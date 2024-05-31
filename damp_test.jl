

#considering the comparison with bowden Em, with and without damping/denisty etc.

using MID
using Plots; plotlyjs()



N = 2000; 
#rgrid = collect(LinRange(0, 1, N));
rgrid = clustered_grid(N, 0.92, 0.98, 0.2)
geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo, δ=-4e-8, dens=axel_dens); 

grids = init_grids(rgrid=rgrid, mstart=1, mcount=4, nstart=-2, ncount=1);
grids = init_grids(rgrid=rgrid, mstart=2, mcount=2, nstart=-2, ncount=1);
#tae_freq = (0.381 / geo.R0)^2


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=(0.39/geo.R0)^2, reconstruct=true);


reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)


tae_ind = 1
display(ω[tae_ind])
plot_potential(r=rgrid, ϕ=ϕ, ind=tae_ind, pmd=grids.pmd, n=1)

display(imag(ω[tae_ind]^2))



#now consider the damping...
#with N=200, δ=-4e-7, R0=10
#fudam = 0.38785232957063126 - 0.001552151780845707im
#mc2 = 0.3928010001211876 - 0.00205577554242077im
#mc4 = 0.3912421249606025 - 0.0018871990469559274im
#tae and damping are both being drastically overestimated.

#with N=1000, δ=-4e-7, R0=10, clustered_grid(N, 0.92, 0.98, 0.2)
#fudam = 0.38786724957293817 - 0.0015539766458243375im #basically no difference.
#mc2 = 0.39266449565197176 - 0.0021706964964117417im
#mc4 = 0.39111378152963566 - 0.0019807359342886677im
#tae and damping are both being drastically overestimated.
#but damping rate predicted by fudam method is in between our method and Axel's results?

#with N=1000, δ=-4e-8, R0=10, clustered_grid(N, 0.92, 0.98, 0.2)
#fudam = 0.38786946860058735 - 0.0011058276694208577im #decreases the damping significanlty.
#mc2 = 0.3926532806264489 - 0.002129619861492239im
#mc4 = 0.3911032476722936 - 0.001937696298444736im

#with N=2000, δ=-4e-8, R0=10, clustered_grid(N, 0.92, 0.98, 0.2)
#fudam = 0.38786946860058735 - 0.0011058276694208577im #decreases the damping significanlty.
#mc2 = 0.3926538800525755 - 0.002130175367374948im
#mc4 = 0.39110378736344276 - 0.0019381909446406754im

((0.38786724957293817 - 0.0015539766458243375im)*0.99)^2

#now with axel dens

#For R0=10, 
#0.3901 from fudam 
#0.3930 with only m=1, 2
#0.3916 with m=0, 1, 2, 3
#difference seems to be more significant.

#For R0=20, 
#0.3949 from fudam 
#0.3963 with only m=1, 2
#0.3960 with m=0, 1, 2, 3



#unifrom dens.

#For R0=10, with uniform dens.
#0.379989 from fudam 
#0.3812 with only m=1, 2
#0.38029 with m=0, 1, 2, 3


#For R0=10, with uniform dens.
#0.38899 from fudam 
#0.3895 with only m=1, 2
#0.3893 with m=0, 1, 2, 3
#seems to be a similar amount of difference for larger R0.

tae_ind = find_ind(ω, 0.398)



display(ω[tae_ind])
display(imag(ω[tae_ind])/real(ω[tae_ind]))

#converging to double somehow???
display(imag(ω[tae_ind]^2)/2)

tae_freq = (ω[tae_ind] / geo.R0)^2