

#testing possible new q-profiles

using MID
using Plots; plotlyjs()
using FastGaussQuadrature


function new_q(r)
    a = 1.05
    b = 0.355

    q = a + b*r^2
    dq = 2 * b * r

    return q, dq
end

Nr = 100;

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=new_q, geo=geo)#, met=no_delta_metric!); 

rgrid = init_fem_grid(N=Nr);
θgrid = init_sm_grid(start=2, count=2);
θgrid = init_fem_grid(N=10)
ζgrid = init_sm_grid(start=-2, count=1);
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid);


ω_cont = continuum(prob=prob, grids=grids);


plot_continuum(ω=ω_cont, grids=grids)


ξr, wgr = gausslegendre(grids.r.gp); #same as python!
ξθ, wgθ = gausslegendre(grids.θ.gp);

#gets the basis 

S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ = MID.hermite_basis(ξr, ξθ);

So, dSro, dSθo, ddSrro, ddSrθo, ddSθθo = MID.hermite_basis(ξr, ξθ);


display(S)
display(So)

println(S[:, :, :] .- So[:, :, :])