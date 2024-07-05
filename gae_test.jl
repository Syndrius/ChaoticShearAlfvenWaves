
#lets see if we can have a look at a gae!
#found one, should probably incorporate this properly once we understand this!

using MID
using Plots; plotlyjs()

#about as big as direct solving can handle
#for this case, fss is much better, 
Nr = 100
Nθ = 10

rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=0, count=3)
θgrid = init_fem_grid(N=Nθ, pf=1)
ζgrid = init_sm_grid(start=0, count=1)

grids = init_grids(rgrid, θgrid, ζgrid)


geo = GeoParamsT(R0=1000) #cylinder??

function gae_q(r)

    #paper says constant, that may be after the Geometry contribution though!

    q = 1/r
    dq = -1.0/r^2
    return q, dq
end

function gae_dens(r)
    κ = -1 #paper just says less than 0?
    p = 1

    return (1-p*r^2)^κ
end

prob = init_problem(q = gae_q, dens=gae_dens, geo=geo)

ω, ϕ = construct_and_solve(grids = grids, prob = prob, full_spectrum=true);


ϕms = mode_structure(ϕ, grids);
#reconstruct_continuum(ω, ϕ, grids)
reconstruct_continuum(ω, ϕms, grids)

#looks like we have found a gae!
#pretty similar results with both methods, fss is significanlty faster though!
gae_ind = find_ind(ω, 0.517)
tae_ind = find_ind(ω, 0.49)
display(ω[tae_ind])


#looks like a very global mode. We need to understand why a turning point in the continuum allows 
#a global mode to form. This will hopefully help us undertsand tae's.
#still need to understand gaps though!

plot_potential(ϕms, grids, gae_ind, 1)

plot_phi_surface(ϕ, grids, gae_ind)


maxm = Int64(grids.θ.N/2 + grids.θ.pf - 1)

if iseven(grids.θ.N)
    mlist = -maxm+1:1:maxm
else
    mlist = -maxm:1:maxm
end
mlist = 0:1:grids.θ.N-1