
#test with island coordinates
using MID
using MIDViz
#%%
geo = init_geometry(:κ, R0=1.0)

isl = init_island(:κ, m0=1, n0=-1, w=0.1, ψ0=0.5, qp=1.0)

#need to change this to allow a single isl to be input.
fields = init_fields(:κ, q=island_q, isl=isl)
prob = init_problem(geometry=geo, fields=fields)
prob.fields.isls[1]
#%%
mt = MID.MetT()
prob.geo.met(mt, 0.1, 0.1, 0.01, 1.0)
mt.dJ
#%%

#should make this grid auto change the start/stop.
κgrid = init_grid(:κ, 20, start=0.001, stop=0.9999)
#κgrid = init_grid(:cont, 30, start=0.001, stop=0.9)
#ᾱgrid = init_grid(:sm, 2, start=0)
ᾱgrid = init_grid(:ᾱ, 5, pf=1)
τgrid = init_grid(:sm, 1, start=0)

grids = init_grids(κgrid, ᾱgrid, τgrid)
#%%
solver = init_solver(prob=prob, full_spectrum=true)
#%%
evals, ϕ, ϕft = compute_spectrum(prob, grids, solver);
#evals = compute_spectrum(prob, grids);
#%%
#cool, shows basic continuum and global mode!
continuum_plot(evals, ylimits=(0.0, 0.05))

ind = find_ind(evals, 0.015)
evals.ω[ind]

harmonic_plot(ϕft, grids, ind)

gap = 0.008


mindiff = 10.0

for ω in evals.ω
    diff = minimum(abs.(gap .- ω))
    if diff < mindiff
        global mindiff = diff
    end
end

mindiff
