

#can be devide a continuum, such that we have an n=2 tae at (say) r=0.5, while having an n0=2 island at another gap at (say) r=0.75?

using MID
using Plots; plotlyjs()
using Printf


#first step is to see if we can create a good continuum!
function test_q(r)

    a = 1.05
    b = 0.8

    q = a+b*r^2
    dq = 2*b*r
    return q, dq
end




Nr = 100;
geo = GeoParamsT(R0=10.0)

prob = init_problem(q=test_q, geo=geo)#, met=no_delta_metric!); 
θgrid = init_sm_grid(start=-5, count=11)
ζgrid = init_sm_grid(start=-5, count = 5, incr=1)
grids = init_grids(Nr, θgrid, ζgrid);

#this does what we want!
ω_cont = continuum(prob, grids);

new_plot_continuum(ω_cont, grids)



ω_cont = continuum(prob = prob, grids=grids);

plot_continuum(ω=ω_cont, grids=grids)




#############


function new_plot_continuum(ω, grids; filename=nothing, ymin=-0.05, ymax=1.05)

    p = scatter(ylimits=(ymin, ymax))
    rgrid, _, mlist, _, _, nlist, _ = MID.instantiate_grids(grids)
    if rgrid[1]==0
        #same condition used when recontsructing
        rgrid = rgrid[2:end]
    end
    rgrid_plot = repeat(rgrid, 1, size(ω)[2])
    for (i, n) in enumerate(nlist)
        scatter!(vcat(rgrid_plot...), vcat(ω[:, :, i]...), label=@sprintf("n=%d", n))
    end

    display(p)

    if !isnothing(filename)
        savefig(p, filename)
    end

end


Nr = 100;

geo = GeoParamsT(R0=10.0)
prob = init_problem(q=island_damping_q, geo=geo);

#start by finding the freq without island.

θgrid = init_sm_grid(start=-10, count=21)
ζgrid = init_sm_grid(start=-8, count = 5, incr=2)
grids = init_grids(Nr, θgrid, ζgrid);

#this does what we want!
ω_cont = continuum(prob, grids);

new_plot_continuum(ω_cont, grids)

rgrid, _, _, _, _, _, _ = MID.instantiate_grids(grids)

size(ω_cont)

scatter(rgrid, ω_cont[:, :, 1], ylimits=(-0.05, 1.05), color=:red)
scatter!(rgrid, ω_cont[:, :, 2], ylimits=(-0.05, 1.05), color=:blue)
scatter!(rgrid, ω_cont[:, :, 3], ylimits=(-0.05, 1.05), color=:green)
scatter!(rgrid, ω_cont[:, :, 4], ylimits=(-0.05, 1.05), color=:orange)
scatter!(rgrid, ω_cont[:, :, 5], ylimits=(-0.05, 1.05), color=:black)

plot_continuum(ω=ω_cont, grids=grids)


""" Just to have!
function island_damping_q(r)
    a = 1.15
    b = 0.4
    q = a + b * r^2
    dq = 2 * b * r
    return q, dq
end
"""