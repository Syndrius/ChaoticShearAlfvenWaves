
#this case has a clearly wrong q-profile or something
#but it gives the difference between using the full equaiton 
#vs using Berks large aspect ration case
#Based on pictures etc, we should be able to replicate this.
#in particular it is supposed to be a rep of Axel's q so we can hopefully work it out.


#still don't trust it!
#total garbage, continuum is wrong, not sure wot the hek this profile is
#also plot of iota start from ~1.05, cleary cannot be this profile...

function calculation_q(r)

    ι = 0.95016-0.67944*r^2+0.62286*r^4-0.41244*r^6+0.1219*r^8+0.0042185*r^10-0.0013979*r^12

    dι = -2*0.67944*r+4*0.62286*r^3-6*0.41244*r^5+8*0.1219*r^7+10*0.0042185*r^9-12*0.0013979*r^11

    q = 1/ι

    dq = -dι/ι^2
    return q, dq
end


using MID
using Plots; plotlyjs()

#this is giving ~1.8 times what we expect.

N = 2000; 
#the collect is a bit annoying, but ok because we will typically use a clustered grid.
#rgrid = collect(LinRange(0, 1, N));

#getting some weird spikes, may need to double check this is working as expected??
rgrid = clustered_grid(N, 0.85, 0.95, 0.25)

geo = GeoParamsT(R0=10.0)

#isl = IslandT(A=4e-5, m0=5, n0=4)

#pretty confident it is convergeing to 0.32780760640103157 - 0.006921689729791451im
#giving ratio as -0.021115097986236567, so consistently above what the literature is giving!
#uses axel_dens in this case.
prob = init_problem(q=calculation_q, geo=geo, δ=-4.0e-9, dens=axel_dens); #probbaly should use geo if it is part of prob,
#prob = init_problem(q=singular_bowden_q, geo=geo, δ=-4e-9, dens=bowden_singular_dens); #probbaly should use geo if it is part 
#even if it is not really used.
grids = init_grids(rgrid=rgrid, mstart=2, mcount=2, nstart=-2, ncount=1);
#tae_freq = (0.381 / geo.R0)^2


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=(0.39/geo.R0)^2, reconstruct=true);


reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)

tae_ind = 1
display(ω[tae_ind])