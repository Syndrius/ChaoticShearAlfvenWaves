
#basis case of computing the continuum

using MID
using Plots; plotlyjs()
using Printf


Nr = 100;
geo = GeoParamsT(R0=10.0)

prob = init_problem(q=island_damping_q, geo=geo)#, met=no_delta_metric!); 

#rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=-12, count=40)
ζgrid = init_sm_grid(start=-12, count=30)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(Nr, θgrid, ζgrid)


ω_cont = continuum(prob, grids);

plot_continuum(ω_cont, grids, ymax=1)

function anal_cont(grids)

    rgrid, _, mlist, _, _, nlist, _ = MID.instantiate_grids(grids)
    #display(collect(rgrid))
    q = zeros(length(rgrid))
    for (i, r) in enumerate(rgrid)
        tq, _ = Axel_q(r)
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

    savefig(p, "label_test.png")

end

anal_cont(grids)