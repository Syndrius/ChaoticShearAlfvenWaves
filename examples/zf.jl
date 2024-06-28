
using MID
using Plots; plotlyjs()
using FFTW

grids = init_grids_zf(Nr=100, Nθ=20, m=2, nstart=-2, ncount=1)

geo = GeoParamsT(R0=4.0)

prob = init_problem(q=Axel_q, geo=geo)#, met=no_delta_metric!);

#construct is significantly longer now, even 200x30 took ~10s.
#but solve is v quick. This may make parallisation more worth, hopefully memory is much smaller.
W, I = construct_zf(prob=prob, grids=grids);

fake_grid = init_grids(N=50, sep1=0.91, sep2=0.98, frac=0.25, mstart=2, mcount=2, nstart=-2, ncount=1);

#ω, efuncs = full_spectrum_solve(Wmat=W, Imat=I, grids=fake_grid, reconstruct=false, R0=geo.R0)

ω, efuncs = arpack_solve(Wmat=W, Imat=I, grids=fake_grid, efuncs=true, σ=(0.37^2)/geo.R0^2, reconstruct=false, R0=geo.R0, nev=50);

#@doc arpack_solve

#actually looks v promising!!!! wild 
#tae freq may be a bit low but need to check first
#need to look at efuncs and reconstruct cont etc.
scatter(ones(length(ω)), real.(ω), ylimits=(-0.05, 1.05))

tae_ind = find_ind(ω, 0.37)

display(ω[tae_ind])

ϕ = reconstruct_phi_zf(efuncs, length(ω), grids);

#looks ok, not particulary periodic though! 
#seems periodic ish, unsure if we maybe need to increase res or perhaps remove last point in θgrid
#or if something is actually wrong!
rgrid = LinRange(0, 1, grids.rd.N)
θgrid = LinRange(0, 2π, grids.θd.N+1)[1:end-1]

#why does the tae surface look so different?
#this must mean our understanding of the solutions is cooked
#maybe we need to multiply by -m in construct_surface??
#tae mode structure in og case does match two mode case,
#but is that because an extra fourier transform was required in that case as well??
surface(θgrid, rgrid, real(ϕ[tae_ind+1, :, :, 1]))

edge1 = real(ϕ[tae_ind, 1, :, 1])
edge2 = real(ϕ[tae_ind, end, :, 1])

#so with Nθ=30, result is periodic up to e-18, not sure if it should be perfect though!
display(maximum(edge1 .- edge2))

ϕms = mode_structure_zf(ϕ, grids)

p = plot()

for i in 1:grids.θd.N
    plot!(rgrid, real.(ϕms[tae_ind, :, i, 1]))
end
display(p)

reconstruct_continuum_zf(ω, ϕms, grids)#, nothing, -0.05, 10)

z = construct_surface()


surface(real(ϕms[tae_ind, :, :, 1]))
#we seem to be getting slightly different results, always good.
#perhaps this is better? Seems a little closer for the fu-dam result. 
#we probably have to paralise for this method to be remotely practical.

#note that this requires that the modes structure already be taken.


    