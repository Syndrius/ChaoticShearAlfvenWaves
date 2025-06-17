
geo = init_geo(R0=4.0)

prob = init_problem(type=:flux, q=MID.Equilibrium.flux_fu_dam_q, geo=geo);

ψgrid = init_grid(type=:rc, N=20, stop=0.5)
θgrid = init_grid(type=:as, start=1, N=2)
ζgrid = init_grid(type=:as, start=-1, N=1)


grids = init_grids(ψgrid, θgrid, ζgrid)


ω = compute_continuum(prob, grids);

#continuum_plot(ω, grids)

gap = 0.32


mindiff = 10.0

for i in 1:θgrid.N, j in 1:ζgrid.N
    diff = minimum(abs.(gap .- ω[:, i, j]))
    if diff < mindiff
        global mindiff = diff
    end
end

#this tests asserst that there is a gap in the continuum in the correct spot.
@test mindiff > 0.04



