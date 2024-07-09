
using MID
using MIDViz
using Plots; plotlyjs()


Nr = 50
Nθ = 15

function slab_island_q(r)
    a = 1
    b = 2
    q = a+b*r
    #not sure about deriv
    dq = b
    return q, dq
end


q = island_damping_q

rgrid = init_fem_grid(N=Nr)

θgrid = init_fem_grid(N=Nθ, pf=2)
#θgrid = init_sm_grid(start=2, count=4)

ζgrid = init_sm_grid(start=-2, count = 1, incr=2)

grids = init_grids(rgrid, θgrid, ζgrid)

isl = IslandT(A=5.0e-4, m0=5, n0=4);

geo = GeoParamsT(R0=100.0)

prob = init_problem(q=q, geo=geo, isl=isl, met=slab_metric!); 
#prob = init_problem(q=slab_island_q, geo=geo, isl=isl, met=slab_metric!); 



ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true, σ = 0.3, nev=100);

length(ω)
#ϕms = mode_structure(ϕ, grids);
omdata = reconstruct_slab_continuum(ω, ϕ, grids)

ind = find_ind(omdata, 0.435)
omdata[ind]

plot_potential(ϕ, grids, ind)

contour_plot(ϕ, grids, ind)

Ntraj = 40
flux_list = LinRange(0.01, 0.5, Ntraj)
rlist = @. sqrt(2*flux_list)
rp, θp = poincare_plot(q, slab_to_plot, 500, Ntraj, 0, 0.2, 4.7, 0, geo.R0, isl, rlist)

plot_contour_poincare(ϕ, grids, ind, rp, θp)