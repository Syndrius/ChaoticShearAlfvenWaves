
geo = init_geometry()
fields = init_fields()


prob = init_problem(fields=fields, geometry=geo)

ψgrid = init_grid(:cont, 20)
θgrid = init_grid(:sm, 2, start=1)
ζgrid = init_grid(:sm, 1, start=-1)


grids = init_grids(ψgrid, θgrid, ζgrid)


evals = compute_spectrum(prob, grids);

#continuum_plot(ω, grids)

gap = 0.32


mindiff = 10.0

for ω in evals.ω
    diff = minimum(abs.(gap .- ω))
    if diff < mindiff
        global mindiff = diff
    end
end

#this tests asserst that there is a gap in the continuum in the correct spot.
@test mindiff > 0.04



