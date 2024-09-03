

using MID
using Plots; plotlyjs()
using Printf


function calculation_q(r)

    ι = 0.95016-0.67944*r^2+0.62286*r^4-0.41244*r^6+0.1219*r^8+0.0042185*r^10-0.0013979*r^12

    dι = -2*0.67944*r+4*0.62286*r^3-6*0.41244*r^5+8*0.1219*r^7+10*0.0042185*r^9-12*0.0013979*r^11

    q = 1/ι

    dq = -dι/ι^2
    return q, dq
end


Nr = 100;
geo = GeoParamsT(R0=10.444)

prob = init_problem(q=calculation_q, geo=geo, dens=axel_dens)#, met=no_delta_metric!); 

#rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=2, count=2)
ζgrid = init_sm_grid(start=-2, count=1)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(Nr, θgrid, ζgrid)


ω_cont = continuum(prob, grids);

plot_continuum(ω_cont.^2, grids, ymax=1)



function anal_cont(grids)

    rgrid, _, mlist, _, _, nlist, _ = MID.instantiate_grids(grids)
    #display(collect(rgrid))
    q = zeros(length(rgrid))
    for (i, r) in enumerate(rgrid)
        tq, _ = default_island_q(r)
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

anal_cont(grids)