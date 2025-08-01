
#using more exact solutions to determine why our interpolation is not working
using MID
using MIDViz
using Plots
#%%

flr = MID.Structures.FLRT(δ=0.1)
prob = MID.Helmholtz.TestProblemT(geo=init_geo(R0=1.0), flr=flr)

rgrid = init_grid(type=:rf, N=20)
θgrid = init_grid(type=:af, N=3)
ζgrid = init_grid(type=:af, N=3)
θgrid = init_grid(type=:as, start=1, N=3)
ζgrid = init_grid(type=:as, start=-2, N=2)

grids = init_grids(rgrid, θgrid, ζgrid)

solver = init_solver(full_spectrum=true, prob=prob)
evals, ϕ, ϕft = compute_spectrum(grids=grids, prob=prob, deriv=true, solver=solver);
#%%

ind1 = argmin(abs.(real.(evals.ω) .- 5.232))
ind2 = argmin(abs.(real.(evals.ω) .- 14.829))

#why out of order? Can only be a bad sign.
#coolio, so it looks like we have some more bugs.
#I thought we go this actually working? guess not.
#more issues to fix, hopefully this is just a scale problemo.
evals.ω #wot the fuck are these eigenvlaues.
evals.ω[ind1]

ind1 = 1963
evals.ω[ind1]
Nmr = 100
Nmθ = 15
Nmζ = 9
mr = LinRange(0, 1, Nmr)
mθ = LinRange(0, 2π, Nmθ)
mζ = LinRange(0, 2π, Nmζ)

rg, θg, ζg = MID.inst_grids(grids)
intphi = zeros(ComplexF64, Nmr)
#intphi2 = zeros(ComplexF64, Nmr)
for i in 1:Nmr
    intphi[i] = MID.Mapping.hermite_interpolation(mr[i], ϕft[ind2, :, 2, 2, :], rg)
    #intphi2[i] = itp(mr[i])
end
#%%
potential_plot(ϕft, grids, ind2)
plot!(mr, real.(intphi))

#shape is good, but the scale is cooked beyond belief.
plot(rg, real.(ϕft[ind2, :, 2, 2, 1]))
plot!(rg, real.(ϕft[ind1, :, 2, 2, 2]))
