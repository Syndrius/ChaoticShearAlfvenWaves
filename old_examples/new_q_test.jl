

using MID
using Plots; plotlyjs() 


#this may be the best option so that the island and everything is at 0.5.
function new_q(r)
    a = 1.15
    b = 0.4
    q = a + b * r^2
    dq = 2 * b * r
    return q, dq
end

N = 100;
#the collect is a bit annoying, but ok because we will typically use a clustered grid.
rgrid = collect(LinRange(0, 1, N));

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=new_q, geo=geo); 
grids = init_grids(rgrid=rgrid, mstart=2, mcount=2, nstart=-2, ncount=1);
#tae_freq = (0.395845)^2/10^2


ω_cont = continuum(prob=prob, grids=grids)

scatter(rgrid[2:end], ω_cont, ylimits=(-0.05, 1.05), legend=false)
plot_continuum(ω = ω_cont, rgrid=rgrid)



ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true);

reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)

tae_ind = find_ind(ω, 0.3828)
plot_potential(r=rgrid, ϕ=ϕ, ind=tae_ind, pmd=grids.pmd, n=1)