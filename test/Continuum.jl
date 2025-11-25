#testing flux case.
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

#testing geometric radius case.
geo = init_geometry()
fields = init_fields(:r)

prob = init_problem(geometry=geo, fields=fields)

rgrid = init_grid(:cont, 20)
θgrid = init_grid(:sm, 2, start=1)
ζgrid = init_grid(:sm, 1, start=-1)


grids = init_grids(rgrid, θgrid, ζgrid)


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


#testing island coordinates
geo = init_geometry(:κ, R0=1.0)

isl = init_island(:κ, m0=1, n0=-1, w=0.1, ψ0=0.5, qp=1.0)

fields = init_fields(:κ, q=island_q, isl=isl)
prob = init_problem(geometry=geo, fields=fields)

#stop the grid early to prevent gap closing at the separatrix.
κgrid = init_grid(:cont, 30, start=0.001, stop=0.9)
ᾱgrid = init_grid(:sm, 2, start=0)
τgrid = init_grid(:sm, 1, start=0)

grids = init_grids(κgrid, ᾱgrid, τgrid)

evals = compute_spectrum(prob, grids)

gap = 0.008


mindiff = 10.0

for ω in evals.ω
    diff = minimum(abs.(gap .- ω))
    if diff < mindiff
        global mindiff = diff
    end
end

#this tests asserst that there is a gap in the continuum in the correct spot.
#note that this gap is smaller!
@test mindiff > 0.005
