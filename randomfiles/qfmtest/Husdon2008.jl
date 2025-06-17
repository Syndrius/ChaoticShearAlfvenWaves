
#trying to replicate the image in Hudson and Breslau 2008.
#this will require modifiying our magnetic field I think!
#could probably get pretty close with the flux version and a slab metric I think!

#this will require a pretty serious shift on how our system is initialised.
#this is probably an overdue change.
#hudson uses iota, naturally...

#think we can replicate the same idea with q, just need to flip the rationality of the islands?

using MID
using MIDViz

#%%

function hudson_q(ψ::Float64)

    #real guess, probably wrong!
    return 1/ψ, -1/ψ^2
end
#%%
#shouldn't matter for slab case.
geo = init_geo(R0=10.0)

k = 2.1e-3
isl1 = init_island(m0=2, n0=-1, A=k/2, flux=true)
isl2 = init_island(m0=3, n0=-2, A=k/3, flux=true)

isls = [isl1, isl2]
#pretty fukin annoying tbh.
#but this doesn work, the above list does not
isls = MID.Geometry.IslandT[]
push!(isls, isl1)
push!(isls, isl2)
isls = MID.Geometry.IslandT[isl1, isl2] #this also works. also shit house. #may need an initialise islands function
#but I think we are getting to the stage where our inputs need to be redone from the ground up... again.
prob = init_problem(met=:slab, B=MID.flux_compute_B!, geo=geo, isls=isls, q=hudson_q)
#%%

Ntraj = 200;
#rlist = collect(LinRange(0.45, 0.65, Ntraj));
rlist = collect(LinRange(0.4, 0.7, Ntraj));
Nlaps = 500;

#ok so pretty sure this matches now,
#unsure why nested parts are not showing?
#key difference is the symmetry.
#unsure how much that matters
#does looks like we have the wrong symmetry or something
#i assume this is just the starting point or the slice we are choosing
#most of the features are still there though!
poincare_plot(prob, Nlaps, Ntraj, rlist; ylimits=(0.5, 0.68))

#%%
#odd we have to flip this
rats1 = lowest_rationals(10, prob.q(0.68)[1], prob.q(0.5)[1])
#will this work given we have a flipped case? probably not?
gl1 = surface_guess(rats1, prob.q)
#looks similar, ours are much smoother. Unsure if this is good or bad,
#probably a good thing for the coordinate transformation.
#changing N to 32 made no difference, our surfaces just look significantly smoother.
#but otherwise they are ok.
#looks like the problem is more the parity? of the islands, the paper has the islands aligned at θ=0.
#this probably just shows us that ours are fine.
surfs1 = construct_surfaces(rats1, gl1, prob, M=32, N=32);
plot_surfs(surfs1)

