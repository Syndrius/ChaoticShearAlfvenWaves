#TODO
#so we need to see if we can figure out why our results are a bit wrong with many m modes
#need to add the derivatives and test the ffs and fff cases
#May need to see if we can add the non-n=0 cases, unsure if this will be practical without a major shift in the code.
#then obvs clean up and check to make sure the coordinate transformation is working etc
#this file is in the wrong place!


using MID
using MIDIslands
using Plots
using MIDViz

#%%

#to match case on gadi!
isl = init_island(m0=1, n0=-1, w=0.1, qp=1.0, ψ0=0.5, coords=true)
#to match zhisongs benchamrk case
#isl = init_island(m0=5, n0=-2, A=1e-5, qp=4.0, ψ0=0.125, coords=true)
isl = MID.inst_island(isl)
isl.A
#%%

#ok, we are computing, probs garbage.
#but this is enough to check if the continuum is computed accuratly
#jks need Jacobian, and maybe its deriv?? -> if we don't need deriv we can get away with taking the determinant
#unsure what the continuum calculation actually requires.
κt = 0.3336
αt = 5.256
τt = 0.0
MID.Geometry.island_metric!(met, κt, αt, τt, 1.0, isl)
#MID.Geometry.old_island_metric!(met, κt, αt, τt, 1.0, isl)

met.gu
#%%
#compare continuums.

#Looks like the problemo is that we are considering the wrong domain
#zeta should be from 0 to 2π*m0. so our fourier expansion is actually different
#unclear how to fix this...
θgrid = init_grid(type=:as, N=15, start=-7)
ζgrid = init_grid(type=:as, N=1, start=0)
#using the below grid with Zhisongs code, matches our case perfectly
#really shows that the domain is the problemo
#unsure how to fix...
#perhaps this is why Axel (and probably us!) only considers n=0?
#for n=0, this does work pretty well
#however, it is still not completly perfect.
#seems to be completly perfect for large κ but not as good for small κ
#may be that our x-axis is a bit wrong
#fixed a typo which made it better, still not perfect like we expect
#ζgrid = init_grid(type=:as, N=1, start=2*isl.m0, incr=isl.m0)

#ok so now it looks like if we just consider a single branch it is absolutely perfect
#but when we do multiple
#they start to veer off, the discrepency is entirely dictated by the starting point
#not by the actual branch
#so not the metric...
#perhaps the alg?
#very small issue, but it would be good to fix.

κlist = LinRange(0.0001, 0.9999, 100)
χlist = @. isl.A - κlist * 2 * isl.A 

ψland = PsiIslandT(isl.m0, isl.n0, isl.A, isl.q0, isl.qp, isl.ψ0)

geo = init_geo(R0=1.0)

quom = island_continuum(χlist, θgrid, ζgrid, geo, ψland, 0)

κp = repeat(κlist, 1, 3)
scatter(κp, sqrt.(abs.(quom .* geo.R0^2)), legend=false, markersize=1.0)#, ylimits=(0, 0.08))

#%%
κc = init_grid(type=:rc, N=100, start=0.0, stop=0.999)
κgrid = init_grid(type=:rf, N=30, start=0.0, stop=0.99999, left_bc=false)
grids = init_grids(κgrid, θgrid, ζgrid)
gridsc = init_grids(κc, θgrid, ζgrid)

prob = init_problem(geo=geo, q=island_q, type=:island, isl=isl, met=:island)

solver = init_solver(full_spectrum=true, prob=prob)

evalsc = compute_continuum(prob, gridsc)

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver)

#looks pretty good, who knows if the global modes are ok though, may need to compare with Axel's case.
continuum_plot(evals, ylimits=(-0.01, 0.03))

scatter!(κp, sqrt.(abs.(quom .* geo.R0^2)), legend=false, markersize=1.0, ylimits=(0, 0.1))

#%%
#seems to be fine now???
#think we still can't do cases without n=0 but that shouldn't matter tbh!
#wot.
#perhaps we didn't reinstantiate the island properly?
κgrid = init_grid(type=:rf, N=40, start=0.0, stop=0.99999, left_bc=false)
#θgrid = init_grid(type=:af, N=10, pf=1)
θgrid = init_grid(type=:as, N=15, start=-7)
ζgrid = init_grid(type=:as, N=1, start=0)
grids = init_grids(κgrid, θgrid, ζgrid)
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver)
continuum_plot(evals, ylimits=(-0.0, 0.1))
hline!([-isl.n0*sqrt(isl.A*isl.qp)])

#%%
#more direct compariosn
harmonic_plot(ϕft, grids, 10)

#%%
κt = 0.3336
αt = 5.256
τt = 0.4
χt = isl.A - κt * 2 * isl.A

α = MIDIslands.Continuum.compute_α(χt, αt, ψland, 0)
ψ = MIDIslands.Continuum.compute_ψ(χt, α, αt, ψland, 0)

fc_met = MID.MetT()

MID.Geometry.flux_cylindrical_metric!(fc_met, ψ, α+τt/isl.q0, τt, geo.R0)

∇χ2 = MIDIslands.Continuum.compute_∇χ2(ψ, α, fc_met, ψland)

MID.Geometry.island_metric!(met, κt, αt, τt, 1.0, isl)
met.gu[1, 1]

#so these are basically identical!
#good news!
∇χ2 / (2*isl.A)^2
#the rest makes use of magnetic field, which doesn't seem to make sense.
#i.e. it is computed from the original B field?
#i guess it only uses the magnitude, which might be the same under a canonical transformation?
#think we might need to compare directly with Zhisong from here on out.
#it is plausible that the minor differences we are seeing now are because our original implementation was wrong

#%%
#trying to read the fkn data from zhisongs code.
using DelimitedFiles

data = readdlm("/Users/matt/phd/MID/convergence/isl_data.csv", ',')


#ok so our MIDIslands is the same as Zhisongs form, meaning, that our metric is a bit wrong!
scatter!(κp, data, legend=false, markersize=1.0, ylimits=(0, 0.4))
