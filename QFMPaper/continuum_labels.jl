
#temp file for getting the branch labels on gadi correct!
using MID
using MIDViz
using Plots; plotlyjs()
#%%

geo = init_geometry(:cyl, R0=1.0)

fields = init_fields(:ψ, q=cantori_q)

prob = init_problem(geometry=geo, fields=fields)
#%%

ψgrid = init_grid(:cont, 200)
θgrid = init_grid(:sm, 10, start=0)
φgrid = init_grid(:sm, 3, start=-6)

grids = init_grids(ψgrid, θgrid, φgrid)
#%%

evals = analytical_spectrum(prob, grids)

continuum_plot(evals, ylimits=(0.11, 0.41), xlimits=(0.5, 0.666), legend=false)
