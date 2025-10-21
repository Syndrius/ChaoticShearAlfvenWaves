
#think we need to get this working again as this would be a good example/test case.
#is it possible/practical to also map the derivatives?
#could be a bit of a disaster?
#think we would need derivatives of all of our mapping equations, which is extremely impractical.
#also not worth the effort for such a small thing!
#but at least getting it to work by mapping the island fake mode to rad/flux demonstrates how it works.
#and should offer some verification that it at least sort of works!
using MID
using MIDViz
using Plots; plotlyjs()
using Plots; gr()
using FFTW
#%%

function global_island_mode(κ, ᾱ, φ)

    m=4
    n=3
    if κ > 1
        return 0
    end
    
    amp = -4*(κ - 0.5)^2 #+1
    amp = sin(π*κ)

    #amp = κ
    dampdκ = -8*(κ-0.5)
    dampdκ = π*cos(π*κ)
    f = amp*exp(1im * m * ᾱ + 1im*n*φ)
    #f = amp*cos(m * ᾱ)
    dfdκ = dampdκ * exp(1im*m*ᾱ + 1im*n*φ)
    dfdᾱ = 1im*m*amp*exp(1im * m * ᾱ + 1im*n*φ)
    dfdφ = 1im*n*amp*exp(1im * m * ᾱ + 1im*n*φ)
    dfdκdᾱ = 1im*m* dampdκ * exp(1im*m*ᾱ + 1im*n*φ)
    dfdκdφ = 1im*n* dampdκ * exp(1im*m*ᾱ + 1im*n*φ)
    dfdᾱdφ = -amp*m*n*exp(1im*m*ᾱ + 1im*n*φ)
    dfdκdᾱdφ = -dampdκ*m*n*exp(1im*m*ᾱ + 1im*n*φ)

    return [f, dfdκ, dfdᾱ, dfdφ, dfdκdᾱ, dfdκdφ, dfdᾱdφ, dfdκdᾱdφ]

end

function cont_island_mode(κ, ᾱ, φ)

    m=2
    if κ > 1
        return 0
    end
    #TODO
    f = exp()*exp(1im * m * ᾱ)
    dfdκ = -8 * (κ - 0.5) * exp(1im*m*ᾱ)
    dfdᾱ = 1im*m*(-4*(κ-0.5)^2+1)*exp(1im * m * ᾱ)
    dfdφ = 0.0
    dfdκdᾱ = -8*1im*m* (κ - 0.5) * exp(1im*m*ᾱ)
    dfdκdφ = 0.0
    dfdᾱdφ = 0.0
    dfdκdᾱdφ = 0.0

    return [f, dfdκ, dfdᾱ, dfdφ, dfdκdᾱ, dfdκdφ, dfdᾱdφ, dfdκdᾱdφ]

end

function create_island_mode(isl_grids)

    κgrid, ᾱgrid, τgrid = MID.inst_grids(isl_grids)
    display(ᾱgrid)

    ϕ_isl = zeros(ComplexF64, length(κgrid), length(ᾱgrid), length(τgrid), 8)

    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, τ) in enumerate(τgrid)
        ϕ_isl[i, j, k, :] .= global_island_mode(κ, ᾱ, τ)
    end
    return ϕ_isl

end
#%%
global_island_mode(0.5, 0.0, 0.0)
ϕf[40, 1, 1, 1]
κi, ᾱi, τi = MID.inst_grids(isl_grids)
κi[40]
ᾱi[1]
τi[1]

contourf(ᾱi, κi, real.(ϕf[:, :, 1, 1]))

size(ϕf)
#%%
ϕplot = MIDViz.Potential.potential_size(ϕf, 1, isl_grids);
z = zeros(ComplexF64, 80, 41)
#z[:, 1:end-1] = ϕplot[:, :, 1] #ϕf[:, :, 1, 1]
#z[:, end] = ϕplot[:, 1, 1] #ϕf[:, 1, 1, 1]
ϕf[argmax(abs.(ϕf[:, :, :, 1])), 1]
z[:, 1:end-1] = ϕf[:, :, 1, 1] / ϕf[argmax(abs.(ϕf[:, :, :, 1])), 1];
z[:, end] = ϕf[:, 1, 1, 1]/ ϕf[argmax(abs.(ϕf[:, :, :, 1])), 1];
x2plot = LinRange(0, 2π, 41)
contourf(x2plot, κi, real.(z))
#%%
#isl = init_island(m0=2, n0=-1, w=0.1, r0=0.5, qp=2.0)
#think this will only ever work properly for (1, -1) island as mapping is not one-to-one
isl = init_island(m0=1, n0=-1, w=0.1, ψ0=0.5, qp=1.0, flux=true)
isl = MID.Geometry.inst_island(isl)
#%%
Nκ = 80
Nᾱ = 40
Nτ = 20
κgrid = init_grid(type=:rf, N=Nκ, stop=0.999)
ᾱgrid = init_grid(type=:af, N=Nᾱ) #probs no pf for now, shouldn't matter
τgrid = init_grid(type=:af, N=Nτ) #probs no pf for now, shouldn't matter
isl_grids = init_grids(κgrid, ᾱgrid, τgrid)

#idealised island mode, created in island geometry.
ϕf = create_island_mode(isl_grids);
contour_plot(ϕf, isl_grids)
ϕf_fft = fft(ϕf, [2, 3]);
harmonic_plot(ϕf_fft, isl_grids)
#%%
#now we transform into toroidal coordinates
Nr = 100 #this might have to be even to prevent the O point from being picked out exactly!
#or an additional guard needs to be put in place.
Nθ = 50
Nζ = 30
rgrid = init_grid(type=:rf, N=Nr, start=0.4, stop=0.6)
θgrid = init_grid(type=:af, N=Nθ)
ζgrid = init_grid(type=:af, N=Nζ)
tor_grids = init_grids(rgrid, θgrid, ζgrid)
#tor_mgrids = MIDIslands.Mapping.MapGridsT(Nx1=Nr, x1max=1.0, Nx2=Nθ, Nx3=Nζ)
#%%
#this is backwards.
#need to reimplement the other one for this to work.
#however, 
#hermite makes this difficult so probabably not worth the effort
#well this doesn't work at all!
isl_to_tor_cm = MID.Mapping.isl_to_tor_coord_map(MID.inst_grids(tor_grids)..., isl);

display(isl_to_tor_cm[8, 1, 1])

ri, θi, ζi = MID.inst_grids(tor_grids)
display((ri[1], θi[1], ζi[1]))

MID.Mapping.tor_coords_to_isl(ri[1], θi[1], ζi[1], isl)
αt = θi[1] - ζi[12] / isl.q0
mod(αt*isl.m0, 2π)
#%%
#should be using ψ now.
ϕtor = zeros(ComplexF64, Nr, Nθ, Nζ); #fkn rip, cannot make this one have the 8. -> wont be able to go the other way. Unless we interpolate the deriv?
#think this might not work because we are not scaling the derivs...
MID.Mapping.efunc_map!(ϕtor, Nr, Nθ, Nζ, ϕf, MID.inst_grids(isl_grids)..., isl_to_tor_cm)
#%%
#try the manual way
#this makes no difference.
itp = interpolate(MID.inst_grids(isl_grids), ϕf[:, :, :, 1], (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));

#seems like this shouldn't be needed!
ext = extrapolate(itp, Periodic());

for i in 1:Nr, j in 1:Nθ, k in 1:Nζ
    #display(tor_to_isl_cm[i, j, k])
    ϕtor[i, j, k] = ext(isl_to_tor_cm[i, j, k]...)
end
#%%


#working now, ideally we could map backwards, but this would require interpolation!
#plausibly could do all the derivs with some wild metric work but by golly that will not be worth the hassle.
ϕtor_fft = fft(ϕtor, [2, 3]);
#this should be symmetric I would think!
#made it way worse lol.
#much happier with the transformation fro isl to tor, but this is still cooked af.
#still weirdly not symmetrical. 
harmonic_plot(ϕtor_fft, tor_grids)
contour_plot(ϕtor, tor_grids)
#%%
#think we can map the otherway, just without hermite obvs.
using Interpolations

Nmκ = 100
Nmᾱ = 60
Nmτ = 30

κmgrid = init_grid(type=:rf, N=Nmκ, stop=0.999)
ᾱmgrid = init_grid(type=:af, N=Nmᾱ) #probs no pf for now, shouldn't matter
τmgrid = init_grid(type=:af, N=Nmτ) #probs no pf for now, shouldn't matter
islm_grids = init_grids(κmgrid, ᾱmgrid, τmgrid)

tor_to_isl_cm = MID.Mapping.tor_to_isl_coord_map(MID.inst_grids(islm_grids)..., isl);
#%%
ϕm_isl = zeros(ComplexF64, Nmκ, Nmᾱ, Nmτ);

itp = interpolate(MID.inst_grids(tor_grids), ϕtor, (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));

#seems like this shouldn't be needed!
ext = extrapolate(itp, Periodic());

for i in 1:Nmκ, j in 1:Nmᾱ, k in 1:Nmτ
    #display(tor_to_isl_cm[i, j, k])
    ϕm_isl[i, j, k] = ext(tor_to_isl_cm[i, j, k]...)
end
#%%

#seems to be working now
#still unsure why the toroidal case is not symmetric???
#think the phase is getting cooked for some reason!
#harmonic structure looks v good though!
#unsure how to fix the phase.
#phase ais a plotting problemo, not a big deal.
ϕm_isl_fft = fft(ϕm_isl, [2, 3]);
harmonic_plot(ϕm_isl_fft, islm_grids)
contour_plot(ϕm_isl, islm_grids)

