
using MID
using Plots; plotlyjs()


#Comparing our code with Axel's contour method,
#In axels case, we significanlty underestimate the damping computed, off by ~43%
#Changing to diagonal metric we are overestimating by ~9%
#increasing R0 to 20
#with off diagonal we get within 2%.

#But consider Bowden singular
#Axel's contour method, with the ϕ equation, does not predict the same result
#rather it predicts -0.018987110268100602, an ~8% difference.
#our code predicts -0.01447212409495936, again underestimating the damping,
#with diagonal metric, we predict -0.018822434541716268, which is almost identical to Axel's contour prediction. and much closer to Bowden Singular.

#Cleary, this buisness is garbage wraped in rubbish
#Seems like everyone's code is rubbish
#big concern for us is that fixing the B error has drastically changed our predictions...
#Seems like our code is good enough though!

#but also why would our code match Contour so well for bowden's but not int he other case???
#perhaps this is because the continuum overlap is closer to the tae peak, meaning the loss of toroidal-ness away from the gap is reduced.





#here we verify the damping calculation by comparing to Bowden and Hole 2015.
#They determine a TAE with normalised frequency of Ω = 0.326 - 0.00571i 
#with a damping ratio of Ωi/Ωr = -0.0175

#this is done by solving the reduced two mode equation originally derived by Berk et al 1992.
#They employ the q-profile q = 1 + 2r^2
#and a density profile of 1/2(1-tanh((r-0.7)/0.05))
#these profiles are designed to close the TAE gap, 
#so that the TAE interacts with the continuum at ~r=0.767, leading to damping.

N = 3000
δ = -4.0e-9
geo = GeoParamsT(R0=10.0)
prob = init_problem(q=singular_bowden_q, geo=geo, δ=δ, dens=bowden_singular_dens)
#the radial grid is clustered around the singularity between 0.75 and 0.8
grids = init_grids(N=N, sep1=0.75, sep2=0.8, frac=0.2, mstart=1, mcount=2, nstart=-1, ncount=1);

tae_freq = (0.326/geo.R0)^2; #un-normalised target frequency.
ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq, reconstruct=true);

#first we verify that we have found a tae.
plot_potential(grids=grids, ϕ=ϕ, ind=1, n=1)

#then print the frequency and damping rate.
display(ω[1])
display(imag(ω[1]) / real(ω[1]))
#The damping rate found is ~22% different than expected.

#This is due to the difference in equation solved.
#Part of Berk et al's derivation neglected the contribution of the off-diagonal terms in the metric.
#if we remove the off diagonal terms from the metric,
prob = init_problem(q=singular_bowden_q, geo=geo, δ=δ, dens=bowden_singular_dens, met=diagonal_toroidal_metric!);
ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq, reconstruct=true);

#verify that we have found a tae.
plot_potential(grids=grids, ϕ=ϕ, ind=1, n=1)

#then print the frequency and damping rate.
display(ω[1])
display(imag(ω[1]) / real(ω[1]))
#this brings us much closer to the expected result, but still a difference of ~7%.


#Another approximation made by Berk et al was to neglect any toroidal corrections
#other than in terms containing second radial derivatives.
#to emulate this, we compute the W and I matrices twice, once using a cylindircal metric
#and the second time using the toroidal metric.
#The final W and I matrices are mostly taken from their cylindrical versions, other than
#the elements that are multiplied by double radial derivatives.
#this is implemented in a separate function, and takes longer due to the extra computations required.
prob = init_problem(q=singular_bowden_q, geo=geo, δ=δ, dens=bowden_singular_dens);
ω, ϕ = analytical_construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq, reconstruct=true);

#verify that we have found a tae.
plot_potential(grids=grids, ϕ=ϕ, ind=1, n=1)

#then print the frequency and damping rate.
display(ω[1])
display(imag(ω[1]) / real(ω[1]))
#now gives a difference of ~3% in damping rate.