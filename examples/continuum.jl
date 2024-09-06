
#basis case of computing the continuum

using MID
using Plots; plotlyjs()
using Printf


Nr = 100;
geo = GeoParamsT(R0=10.0)

prob = init_problem(q=test_q, geo=geo)#, met=no_delta_metric!); 

#rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=2, count=2)
ζgrid = init_sm_grid(start=-2, count=1)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(Nr, θgrid, ζgrid)


ω_cont = continuum(prob, grids);

plot_continuum(ω_cont, grids, ymax=1)


function test_q(r)

    a = 1.2
    b = 8/15
    q = a + b*r^2
    dq = 2 * b * r
    return q, dq
end

display(island_damping_q(1))

function anal_cont(grids)

    rgrid, _, mlist, _, _, nlist, _ = MID.instantiate_grids(grids)
    #display(collect(rgrid))
    q = zeros(length(rgrid))
    for (i, r) in enumerate(rgrid)
        #tq, _ = default_island_q(r)
        tq, _ = Axel_q(r)
        #tq, _ = island_damping_q(r)
        #tq, _ = test_q(r)
        #tq, _ = island_3_2_q(r)
        q[i] = tq
    end
    p = scatter(ylimits=(-0.05, 1))
    #res = Axel_q.(collect(rgrid))
    #q = res[:, 1]
    #display(q)
    for m in mlist
        #display(m)

        for n in nlist
            
            kpar = @. m/q + n
            #display(kpar)
            scatter!(rgrid, abs.(kpar), label=@sprintf("(%d, %d)", m, n))
        end
    end
    display(p)

    #savefig(p, "label_test.png")

end

rid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=-10, count=21)
ζgrid = init_sm_grid(start=-6, count=13)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(Nr, θgrid, ζgrid)

anal_cont(grids)


