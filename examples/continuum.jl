using MID
using Plots
#%%
#maybe give this a default R0 value of 4.0?
#perhaps kwards are useful for these intialisation type things
#unless they are required eg in init_grid, the symbol and N.
geo = init_geometry(4.0)

fields = init_fields()

prob = init_problem(geo=geo, fields=fields)
#%%
rgrid = init_grid(:cont, 20)
θgrid = init_grid(:sm, 2, start=1)
ζgrid = init_grid(:sm, 1, start=-1)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%

evals = compute_spectrum(prob, grids)

scatter(evals.x1, real.(evals.ω))
#%%
gap = 0.32


mindiff = 10.0

for ω in evals.ω
    diff = minimum(abs.(gap .- ω))
    if diff < mindiff
        global mindiff = diff
    end
end

#looks like this test is working a bit.
display(mindiff)

