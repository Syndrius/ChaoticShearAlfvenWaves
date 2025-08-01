
#first actual convergence test with the real case.
#don't think we can realisticly do this, as higher grids results in more
#k vals, so I think we are better off just looking at a tae frequency as the grid resolution goes up.
using MID
using Plots; plotlyjs()
using MIDViz
#%%

flr = MID.Structures.FLRT(δ=1e-8)
geo = init_geo(R0=4.0)
prob = MID.Helmholtz.TestProblemT(flr=flr, geo=geo)

prob = init_problem(type=:flux, q=fu_dam_q, geo=geo)

Nx = 10
Ny = 5
Nz = 3
xgrid = init_grid(type=:rf, N = Nx)#, gp=3)
#ygrid = init_grid(type=:as, start=1, N=3)
#zgrid = init_grid(type=:as, start=-1, N=1)
ygrid = init_grid(type=:af, N=Ny)
zgrid = init_grid(type=:af, N=Nz)

grids = init_grids(xgrid, ygrid, zgrid)

solver = init_solver(full_spectrum=true, prob=prob)
#%%
evals, ϕ, ϕft = MID.Solve.compute_spectrum(prob=prob, grids=grids, solver=solver, deriv=true);
#%%

continuum_plot(evals)

ind1 = find_ind(evals, 0.206)
ind2 = find_ind(evals, 0.1789)

potential_plot(ϕ, grids, ind1, legend=false)
potential_plot(ϕft, grids, ind2)

#%%
Nmx = 50
Nmy = 20
Nmz = 10
mr = init_grid(type=:rf, N=Nmx)
mθ = init_grid(type=:af, N=Nmy)
mζ = init_grid(type=:af, N=Nmz)

mgrids = init_grids(mr, mθ, mζ)

intϕ1 = zeros(ComplexF64, Nmx, Nmy, Nmz);
intϕ2 = zeros(ComplexF64, Nmx, Nmy, Nmz);

xg, yg, zg = MID.inst_grids(grids)
mrg, mθg, mζg = MID.inst_grids(mgrids)
#%%
for i in 1:Nmx, j in 1:Nmy, k in 1:Nmz
    intϕ1[i, j, k] = MID.Mapping.hermite_interpolation(mrg[i], mθg[j], mζg[k], ϕ[41, :, :, :, :], xg, yg, zg)
    intϕ2[i, j, k] = MID.Mapping.hermite_interpolation(mrg[i], mθg[j], mζg[k], ϕ[45, :, :, :, :], xg, yg, zg)
end

ϕ[41, :, :, :, :];
MID.Mapping.hermite_interpolation(mrg[4], mθg[3], mζg[2], ϕ[12, :, :, :, :], xg, yg, zg);
intϕ1
#%%
intϕ1fft = fft(intϕ1, [2, 3]);
intϕ2fft = fft(intϕ2, [2, 3]);
#%%
#wow this is completly fucked.
#how is this possibly so bad. what are we even doing.
p = plot()
for i in 1:size(intϕ1fft)[2], j in 1:size(intϕ1fft)[3]
    plot!(mrg, real.(intϕ1[:, i, j]))
end
display(p)
#%%
#interpolation looks worse, and the fourier transform doesn't work
#so I guess the mapping function probably doesn't work.
contourf(mθg, mrg, real.(intϕ1[:, :, 1]))
contourf(yg, xg, real.(ϕ[41, :, :, 1, 1]))
#%%
#################################################
nbcs = length(unique(MID.Indexing.compute_boundary_inds(grids)))
ktot, code = anal_evals(1:20, ygrid.start:ygrid.start+ygrid.N-1, zgrid.start:zgrid.start+zgrid.N-1);
ktrue = sort(sqrt.(ktot))
diff1 = real.(evals.ω)[nbcs+1:nbcs+20] .- ktrue[1:20]

#%%
real.(evals.ω)[nbcs+1:nbcs+10]
evals.modelabs[nbcs+1:nbcs+10]
ktrue[1:10]
code[3]
diff1[13]
code[7]
collect(zgrid.start:zgrid.start+zgrid.N)

#so actually increasing the number of modes made them more accurate, I guess
#not sure if this is actually true, it may have just meant more low evals are valid
#hence we are only able to find the lowest ones?
#we only really care about fff, so we probably need to start putting this in parallel.
#then just run it with multiple N's
#and probably also check how the phase factor works, probably easiest to label the modes with the mode number
#and hoepflly we see that evals with m=pf converge faster than others?
#I think we need to pick certain evals to focus on, rather than doing a full spread
#doing them all 
scatter(diff1)
scatter!(diff3)
scatter!(diff51)

real.(evals.ω)[nbcs+1:nbcs+20]

