
#changing the og island q-profile to have less islands for better mode resolution

using MID


function new_q(r)
    #so 3/2 and 4 does not work??? unsure why...
    q0 = 3/2
    qp = 1.6
    ψ0 = 0.125
    q = 1 / (1 / q0 - qp / (q0)^2 * (r^2/2-ψ0))
    dq = 4*(q0)^2* qp * r / (2*(q0)-qp*r^2 + 2*qp*ψ0)^2
    return q, dq
end

Nr = 100;
geo = GeoParamsT(R0=5.0)

prob = init_problem(q=new_q, geo=geo)#, met=no_delta_metric!); 

#rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=0, count=14)
ζgrid = init_sm_grid(start=-8, count=8)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(Nr, θgrid, ζgrid)


ω_cont = continuum(prob, grids);

plot_continuum(ω_cont, grids, ymax=1)



using MIDViz


isl = IslandT(m0=3, n0=-2, A=4e-3)


Ntraj = 40
flux_list = LinRange(0.01, 0.5, Ntraj)
rlist = @. sqrt(2*flux_list)
rp, θp = poincare_plot(prob.q, slab_to_plot, 500, Ntraj, 0, 0.0, 0.0, 0.0, prob.geo.R0, isl, rlist);