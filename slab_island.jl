
using MID
using MIDViz
using Plots; plotlyjs()
using LaTeXStrings


Nr = 60
Nθ = 10


q = island_damping_q

rgrid = init_fem_grid(N=Nr)

#pf is defs cooked :(
θgrid = init_fem_grid(N=Nθ)#, pf=2)
#θgrid = init_sm_grid(start=2, count=4)

ζgrid = init_sm_grid(start=-2, count = 1, incr=2)

grids = init_grids(rgrid, θgrid, ζgrid)

isl = IslandT(A=4.0e-4, m0=5, n0=-4);

geo = GeoParamsT(R0=1000.0)

prob = init_problem(q=q, geo=geo, isl=isl)#, met=slab_metric!); 
#prob = init_problem(q=slab_island_q, geo=geo, isl=isl, met=slab_metric!); 


#hmm probably can't see any island modes with only ζ=1...
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true, target_freq = 0.3, nev=100);


#omdata = plot_slab_continuum(evals);
#so ffs matches fss. But it looks like complete garbage in either case... Also slab does not match cylinder...
plot_continuum(evals);

ind = argmin(abs.(omdata .- 0.818))
ind = find_ind(evals, 0.6275764)
omdata[ind]
ind = 2

#so the island seems to be cooking the labelling/mode structure or something. I am v confused...
plot_potential(ϕft, grids, ind)

#why tf is there 6??? surely should be five???? 
#wot the fk is this.
contour_plot(ϕ, grids, ind=ind)

#rgrid, θgrid, _, _, _ = instantiate_grids(grids);

#wot whi is it stuck at r=0??
#contourf(θgrid, rgrid, real.(ϕ[ind, :, :, 1]), levels=100, color=:turbo)



Ntraj = 40
flux_list = LinRange(0.01, 0.5, Ntraj)
rlist = @. sqrt(2*flux_list)
rp, θp = poincare_plot(q, slab_to_plot, 500, Ntraj, 0, 0.2, 4.7, 0, geo.R0, isl, rlist)

plot_contour_poincare(rp, θp, ϕ, grids, ind)

function plot_slab_continuum(evals::MID.EvalsT; filename=nothing, n=nothing, ymin=-0.05, ymax=1.05)
    omdata = real.(evals.ω) ./ evals.r 
    p = scatter(evals.r, omdata, group=evals.modelabs, xlabel=L"r", ylabel=L"\frac{\omega  R_0}{v_A}", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10, ylimits=(ymin, ymax))

    display(p)
    if !isnothing(filename)
        savefig(p, filename)
    end

    return omdata
end