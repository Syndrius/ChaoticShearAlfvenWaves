
using MID




function contrived_q(r)
    m = 4
    n = -1
    r0 = 0.5
    a = 1
    b = 0.2
    #c = (-1-2*m-2*a*n - 2*b*n*r0)/(2*n*r0^2)
    #q = a + b*r * c*r^2
    q = a + 1/(b*r^2+r^4+r^6)
    dq = b+2 #* c * r
    return q, dq
end



Nr = 500;
geo = GeoParamsT(R0=3.0)

prob = init_problem(q=contrived_q, geo=geo)#, met=no_delta_metric!); 

#rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=4, count=2)
ζgrid = init_sm_grid(start=-1, count=1)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(Nr, θgrid, ζgrid)


ω_cont = continuum(prob, grids);

plot_continuum(ω_cont, grids, ymax=2)


rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=4, count=2)
ζgrid = init_sm_grid(start=-1, count=1)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)
evals, ϕ, ϕft = compute_spectrum(grids=grids, prob=prob);

plot_continuum(evals)

tae_ind = find_ind(evals, 0.07)

plot_potential(ϕft, grids, tae_ind)