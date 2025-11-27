
#probably an example of a benchmark qfm use case for a non-resonant island
#then I think also an example of the chaotic case!

#TODO
#non-resonant qfm case
#can also show off mapping here.
#may also want to do some mapping with an island case.

using MID
using MIDViz
#%%

surf_dir = abspath(joinpath(pathof(MID), "../../test/data/"))
surfs = surfaces_from_file(joinpath(surf_dir, "benchmark_surfaces.jld2"));
#%%

geo = init_geometry(:tor, R0=4.0)

#probably cant just use w huh.
#this would require the root solve shite.
isl = init_island(m0=3, n0=1, A=0.1)

fields = init_fields(:ψ, q=cantori_q, isl=isl)

prob = init_problem(geometry=geo, fields=fields)

#%%

sgrid = init_grid(:s, 30, start=0.15, stop=0.9)
ϑgrid = init_grid(:ϑ, 5, pf=1)
ζgrid = init_grid(:ζ, 1, pf=-1)

grids = init_grids(sgrid, ϑgrid, ζgrid)
#%%

solver = init_solver(prob=prob, full_spectrum=true)
#%%

evals, ϕ, ϕft = compute_spectrum(prob, grids, solver, surfs, deriv=true);


continuum_plot(evals)

tae_ind = find_ind(evals, 0.25) 
cont_ind = find_ind(evals, 0.20)

#looks good!
harmonic_plot(ϕft, grids, tae_ind)
harmonic_plot(ϕft, grids, cont_ind)
#%%
#first create grids to map to, 
#note that this can be slow af.

ψgrid = init_grid(:ψ, 80, start=0.25, stop=0.8)
θgrid = init_grid(:θ, 30)
φgrid = init_grid(:φ, 10)
tor_grids = init_grids(ψgrid, θgrid, φgrid)
#%%
#pretty sure this only works for fff.
tor_evals, ϕ_tor, ϕft_tor = qfm_spectrum_to_tor(evals, ϕ, ϕft, grids, tor_grids, surfs);
size(ϕ)

continuum_plot(tor_evals)
tor_ind = find_ind(tor_evals, 0.25)
display(ϕft_tor)
#think these harmonic plots are fkn stoopid.
harmonic_plot(ϕft_tor, tor_grids, tae_ind, label_max=0.5)
harmonic_plot(ϕft_tor, tor_grids, cont_ind, label_max=0.5)
#this is kind of cool
#maybe shows us that our qfm choice was a bit cooked.
#i.e. the perturbation was perhaps a bit large.
contour_plot(ϕ, grids, cont_ind)
contour_plot(ϕ_tor, tor_grids, cont_ind)
