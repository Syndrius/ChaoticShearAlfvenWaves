#maybe this is good enough now
#who the fek knows
#unbelievable that any of this was working
#the mapping is cooked, don't know how the contour plots where
#even kind of working
#shows why the harmonics structrue was so cooked
#don't really understand why our fake cases where working so well
#add this to the pile of stuff.
using MID
using MIDViz
using Plots; plotlyjs()
#%%

prob = init_problem(type=:flux, geo=init_geo(R0=4.0), q=fu_dam_q)

rgrid = init_grid(type=:rf, N=50)
θgrid = init_grid(type=:af, N=5)
ζgrid = init_grid(type=:af, N=3)
θgrid = init_grid(type=:as, start=1, N=3)
ζgrid = init_grid(type=:as, start=-2, N=2)

grids = init_grids(rgrid, θgrid, ζgrid)

solver = init_solver(full_spectrum=true, prob=prob)
evals, ϕ, ϕft = compute_spectrum(grids=grids, prob=prob, deriv=true, solver=solver);
#%%
#small test case, but does show that the interpolation of two efuncs at the same radial location 
#can be split by the interpolation
#this may allow a smoother continuum.
continuum_plot(evals)
ind1 = find_ind(evals, 0.16396)
ind2 = find_ind(evals, 0.1586)

ϕ1 = ϕ[ind1, :, :, :, :];
ϕ2 = ϕ[ind2, :, :, :, :];
ϕ1fft = ϕft[ind1, :, :, :, :];
ϕ2fft = ϕft[ind2, :, :, :, :];

size(ϕ1)
#%%
Nmr = 100
Nmθ = 15
Nmζ = 9
mr = LinRange(0, 1, Nmr)
mθ = LinRange(0, 2π, Nmθ)
mζ = LinRange(0, 2π, Nmζ)

rg, θg, ζg = MID.inst_grids(grids)
intphi = zeros(Nmr, Nmθ, Nmζ);
for i in 1:Nmr, j in 1:Nmθ, k in 1:Nmζ
    intphi[i, j, k] = MID.Mapping.hermite_interpolation(mr[i], mθ[j], mζ[k], ϕ1, rg, θg, ζg)
end
#%%
#given this works for test case in 1d, means it cannot be the order of h's that is wrong, rather
#the actual derivative is wrong? or doesn't match our normal case..
rg2 = range(0, 1, rgrid.N)
itp = cubic_spline_interpolation(rg2, ϕ[ind1, :, 1, 2, 1])
intphi = zeros(ComplexF64, Nmr)
intphi2 = zeros(ComplexF64, Nmr)
for i in 1:Nmr
    intphi[i] = MID.Mapping.hermite_interpolation(mr[i], ϕ[ind1, :, 1, 2, :], rg)
    intphi2[i] = MID.Mapping.hermite_interpolation(mr[i], ϕ[ind2, :, 1, 2, :], rg)
    #intphi2[i] = itp(mr[i])
end

#%%
plot(mr, real.(intphi))
plot!(mr, real.(intphi2))

plot(rg, real.(ϕft[ind1, :, 1, 2, 1]))
#%%
potential_plot(ϕ1, grids)
potential_plot(ϕ2fft, grids)
#%%
p = plot()
for i in 1:Nmθ, j in 1:Nmζ
    plot!(mr, real.(intphi[:, i, j]))
end
display(p)
#%%
MID.Indexing.index_to_grid(2, grids)
MID.Indexing.grid_to_index(1, 1, 1, 1, 1, 2, grids)
#%%
intphifft = fft(intphi, [2, 3]);


p = plot()
for i in 1:Nmθ, j in 1:Nmζ
    plot!(mr, real.(intphifft[:, i, j]))
end
display(p)
#%%
itp = interpolate((rg, θg, ζg), ϕ1[:, :, :, 1], (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));
ext = extrapolate(itp, Periodic());
intphi = zeros(Nmr, Nmθ, Nmζ);
#this is unbeleivably different to the hermite case
#this actually makes the hermite case look good
#obvs not confident it is actually working at all.
for i in 1:Nmr, j in 1:Nmθ, k in 1:Nmζ
    intphi[i, j, k] = ext(mr[i], mθ[j], mζ[k])
end
#%%
basis_id = MID.Indexing.basis_id
grid_id = MID.Indexing.grid_id

h1 = 1
h2 = 1
h3 = 4
bi = 1 + basis_id[h3] + 2 * basis_id[h2] + 4 * basis_id[h1]
