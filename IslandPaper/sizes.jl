#we want to check that the poincare plot is alright!
#and that the gap is relatively large, I think we want to have a small and large 2,-1 case
#and a small and large 3, -2 case.
#we want to consider multiple island case/sizes
#here we will figure out which would be best
#think we will stick with the (2, -1) island
using MID
using MIDViz
using MIDIslands
using Plots; plotlyjs()
using Plots; gr()
function convert_isl(isl::MID.IslandT)
    #converts a mid island (in terms of r)
    #to zhisongs island, in terms of ψ
    #looks to be working!
    qp = isl.qp / isl.r0
    A = (isl.w/4)^2 * qp / isl.q0^2
    ψ0 = isl.r0^2 / 2
    return PsiIslandT(isl.m0, isl.n0, A, isl.q0, qp, ψ0)
end
#%%
#all 4 cases match the contiinuum perfectly! v nice.
#think we can be pretty happy with these 4 island choices

#think this is pretty good, corresponds to ~10% and ~20% of the domain.
isl21a = init_island(m0=2, n0=-1, w=0.05, r0=0.5, qp=2.0)
isl21b = init_island(m0=2, n0=-1, w=0.1, r0=0.5, qp=2.0)

#these widths are the same
#this has a significantly larger gap. Probably a good thing!
#also the gap scales with width. Think there is a relationship somewhere.
isl32a = init_island(m0=3, n0=-2, w=0.05, r0=0.5, qp=2.0) #may need to change qp
isl32b = init_island(m0=3, n0=-2, w=0.1, r0=0.5, qp=2.0) #may need to change qp

isl = isl32b

geo = init_geo(R0=1000.0)
isl_prob = init_problem(q=MID.island_coords_q, met=:island, geo=geo, isl=isl)
tor_prob = init_problem(q=MID.Equilibrium.island_equiv_q, met=:cylinder, geo=geo, isl=isl)
#%%
#poincare plots.

Ntraj = 80
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

x, z = poincare_plot(tor_prob, Nlaps, Ntraj, rlist)

#%%
κgrid = init_grid(type=:rc, N=80, start=0.0, stop=0.999)
ᾱgrid = init_grid(type=:as, N=6, start=-2)
τgrid = init_grid(type=:as, N=3, start=-1)
grids = init_grids(κgrid, ᾱgrid, τgrid)
#%%
ωcont = compute_continuum(isl_prob, grids);
#%%
continuum_plot(ωcont, grids, ylimits=(0.0, 0.3))

#%%
ψland = convert_isl(isl)
pmd = init_grid(type=:as, start=-2, N=6, incr=1)
tmd = init_grid(type=:as, start=0, N=1, incr=2)

#A = isl.A #almost certainly wrong
A = ψland.A #almost certainly wrong
#@. doesn't work with exp of linrange for some reason!
χlist1 =  1 .- 0.01 * exp.( -1 .* collect(LinRange(0, 12, 11)));
χlist2 = LinRange(0, sqrt(χlist1[1]), 191)[1:end-1] .^2 .+ 1e-5;
#ok this matches Zhisong, again no fkn idea what the hell this list is.
χlist = A .- vcat(χlist2, χlist1) .* 2 .* A;
#%%
#may actually be different!
#not ideal, will need to test with less modes. takes like ~5mins.
#ω2list = trapped_continuum(χlist, pmd, tmd, geo, isl);
ω2list = island_continuum(χlist, pmd, tmd, geo, ψland, 0);
#%%
#now that we have a better understanding of the coords, we might want to formalises the transformation between 
#out stuff and Zhisongs, I think any paper should actually use Zhisongs work.
width = 4 * sqrt(ψland.A * ψland.q0^2/isl.qp)
ψ_isl = 2 * width / (π * ψland.m0)
ψ̄m = ψ_isl

ψ̄list = compute_ψ̄(ψland, χlist, 0);
r = repeat(sqrt.(2 .* ψ̄list), 1,  pmd.N * tmd.N);
#%%

#looks to match, x-axis scaling is cooked af.
#as expected!
scatter!(r, sqrt.(abs.(ω2list .* geo.R0^2)), legend=false, markersize=0.1, ylimits=(0, 0.3))
