
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
#%%

function global_island_mode(κ, ᾱ, φ)

    m=2
    if κ > 1
        return 0
    end
    f = (-4*(κ-0.5)^2+1)*exp(1im * m * ᾱ)
    dfdκ = -8 * (κ - 0.5) * exp(1im*m*ᾱ)
    dfdᾱ = 1im*m*(-4*(κ-0.5)^2+1)*exp(1im * m * ᾱ)
    dfdφ = 0.0
    dfdκdᾱ = -8*1im*m* (κ - 0.5) * exp(1im*m*ᾱ)
    dfdκdφ = 0.0
    dfdᾱdφ = 0.0
    dfdκdᾱdφ = 0.0

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

    ϕ_isl = zeros(ComplexF64, length(κgrid), length(ᾱgrid), length(τgrid), 8)

    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, τ) in enumerate(τgrid)
        ϕ_isl[i, j, k, :] .= global_island_mode(κ, ᾱ, τ)
    end
    return ϕ_isl

end
#%%

isl = init_island(m0=2, n0=-1, w=0.1, r0=0.5, qp=2.0)
isl = init_island(m0=2, n0=-1, w=0.1, ψ0=0.5, qp=2.0, flux=true)
isl = MID.Geometry.inst_island(isl)
#%%
κgrid = init_grid(type=:rf, N=Nκ, stop=0.999)
ᾱgrid = init_grid(type=:af, N=Nᾱ) #probs no pf for now, shouldn't matter
φgrid = init_grid(type=:af, N=Nφ) #probs no pf for now, shouldn't matter
isl_grids = init_grids(κgrid, ᾱgrid, φgrid)

#idealised island mode, created in island geometry.
ϕf = create_island_mode(isl_grids);
contour_plot(ϕf, isl_grids)
#%%

#now we transform into toroidal coordinates
Nr = 500
Nθ = 100
Nζ = 10
rgrid = init_grid(type=:rf, N=Nr)
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
tor_to_isl_cm = MID.Mapping.isl_to_tor_coord_map(κf, ᾱf, φf, isl);
#%%
#should be using ψ now.
ϕtor = zeros(ComplexF64, Nr, Nθ, Nζ); #fkn rip, cannot make this one have the 8. -> wont be able to go the other way. Unless we interpolate the deriv?
MID.Mapping.efunc_map!(ϕtor, Nr, Nθ, Nζ, ϕf, κf, ᾱf, φf, tor_to_isl_cm)
#%%


#working now, ideally we could map backwards, but this would require interpolation!
#plausibly could do all the derivs with some wild metric work but by golly that will not be worth the hassle.
potential_plot(ϕtor, tor_grids)
MIDViz.Plotting.contour_plot(ϕtor, tor_grids)
