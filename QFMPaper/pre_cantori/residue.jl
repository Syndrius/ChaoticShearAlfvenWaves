#see how the residue calculations compare with the qfm saw construction.
using MID
using MIDCantori
using MIDViz
using Plots; plotlyjs()
using JLD2
#%%
k_min = 0.0012 #ideally we could get this theoretically!

k_mid = 0.0016

k_max = 0.0019 #above this the qfm surfaces inside the chaotic region get cooked, this may still be too much though!

k = k_min

geo = init_geo(R0=4.0)
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=5, n0=-3, A=k/5, flux=true)

isls = MID.Geometry.IslandT[isl1, isl2]

prob = init_problem(geo=geo, q=MID.Equilibrium.cyl_qfm_q, isls=isls, met=:cylinder, type=:flux)

surfs = load_object("/Users/matt/phd/MID/data/surfaces/flux_qfm_min_surfs.jld2");
#%%
sgrid = init_grid(type=:rf, N = 30, start=0.05, stop=0.95, sep1=0.5, sep2=0.66, frac=0.5)
ϑgrid = init_grid(type=:af, N = 6, pf=2) 
φgrid = init_grid(type=:af, N = 1, pf=-1)
grids = init_grids(sgrid, ϑgrid, φgrid)
#%%
solver = init_solver(nev=150, targets=[0.3], prob=prob)
#%%
evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=surfs, deriv=true);
#%%

#we may have to more carefully choose our s grid, so that we actually pick out the rationals exactly!
continuum_plot(evals, legend=false, xlimits=(0.5, 0.666666))#, n=-2)

ind = find_ind(evals, 0.3571255)

potential_plot(ϕft, grids, ind)

#%%

#we are probbaly going to have to redefien our surfaces, for when most residues are below 0.25
#and for when most are above.
k = 0.0008

geo = init_geo(R0=4.0)
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=5, n0=-3, A=k/5, flux=true)

isls = MID.Geometry.IslandT[isl1, isl2]

prob = init_problem(geo=geo, q=MID.Equilibrium.cyl_qfm_q, isls=isls, met=:cylinder, type=:flux)
#%%
using Plots; gr()
Ntraj = 300;
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
ψlist = collect(LinRange(0.4, 0.8, Ntraj));
poincare_plot(prob, 500, ψlist, zeros(Ntraj), ylimits=(0.5, 0.6666))
#%%
alist, blist = farey_tree(5, 3, 2, 5, 3)

Nrats = length(alist);
Rlist = zeros(Nrats);

for i in 1:Nrats
    #don't love that this one is b/a as the guess. Pretty annoying tbh!
    #think we need to change this to the 1d version.
    #i.e. make it much more specific to the problem we are actually solving!
    guess = surface_guess([(alist[i], blist[i])], prob.q)
    r0, θ0 = MIDCantori.QFM.construct_surfaces([(alist[i], blist[i])], guess, prob, N=4)
    #so this is wrong, 
    _, R = MIDCantori.Residue.residue(r0[1], θ0[1], alist[i], blist[i], prob)
    

    Rlist[i] = R
end
#%%
using Plots; plotlyjs()
guess = surface_guess([(alist[1], blist[1])], prob.q)
xvals = blist ./ alist

#hmmm this is wrong af.
scatter(xvals, clamp.(abs.(Rlist), 0.0, 10.0))



