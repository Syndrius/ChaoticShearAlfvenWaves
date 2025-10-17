
using MIDCantori
using MID
using MIDViz
#%%


k = 0.0015

rationals = [(4, 3), (3, 2), (7, 5), (5, 3), (5, 4), (7, 4), (6, 5), (8, 5), (9, 5), (7, 6), (11, 6), (8, 7), (9, 7), (11, 8), (10, 7)]

gl = 0.6 .* ones(length(rationals))

geo = init_geo(R0=1.0)
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)

isls = MID.Geometry.IslandT[isl1, isl2]

prob = init_problem(geo=geo, q=cantori_q, isls=isls, met=:cylinder, type=:flux)
surfs = construct_surfaces(rationals, gl, prob);


#%%

plot_surfs(surfs)

#%%

ψlist = LinRange(0.5, 0.67, 100)
θlist = zeros(length(ψlist))

poincare_plot(prob, 500, ψlist, θlist, ylimits=(0.5, 0.667))

plot_surfs(surfs, overlay=true, linewidth=2.0, color=:red)

#%%

ψlist = LinRange(0.5, 0.7, 1)
θlist = zeros(length(ψlist))

poincare_plot(prob, 200, ψlist, θlist, surfs=surfs)
#%%
CT = MID.CoordTransformT()

met = MID.MetT()
B = MID.BFieldT()

sitp, sd = MID.create_surf_itp(surfs);

MID.coord_transform!(0.6159, 1.2529, 1.7775, CT, sitp, sd)

display(CT.coords)
