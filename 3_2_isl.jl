
#determining the q-profile to use for a 3/2 island with the usual island q-profile
using MID


function test_island_q(r::Float64)

    q0 = 3.0/2.0
    qp = 0.5
    r0 = 0.5
    #need to verify this!
    #maybe this will fix the flux surface problemo.
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (r^2-r0^2))
    dq = 4*qp*q0^2*r*r0 / (2*q0*r0 - qp * (r^2-r0^2))^2
    return q, dq


end


geo = GeoParamsT(R0 = 1000)

prob = init_problem(q=test_island_q, geo=geo)#, isl=isl)

rgrid = MID.ContGridDataT(N=100)
θgrid = asm_grid(start=0, N=8)
ζgrid = asm_grid(start=-4, N=5)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)


evals = cyl_cont(prob, grids);






#evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);#, target_freq=10);

continuum_plot(evals)#, ymax=0.08)#, ymax