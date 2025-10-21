
#fk the island functions domain basically
using MID
using Plots; gr()
using Plots; plotlyjs()
using Elliptic
#%%

isl = init_island(m0=1, n0=-1, w=0.2, ψ0=0.5, qp=1.0, flux=true)
#wowee the coordinate map seems to be perfect for the radial case???
#how in the damn hek.
#isl = init_island(m0=3, n0=-2, w=0.2, r0=0.5, qp=1.0)#, flux=true)
isl = MID.Geometry.inst_island(isl)
isl.q0
Nκ = 200
Nᾱ = 200

#ok we can actually see evidence of this being completly cooked now!
#still some brokem symmetry somewhere
ζval = 6*π-0.3#3.455751 #3.455751 shows that this is wrong af.

κlist = LinRange(0, 1, Nκ)
ᾱlist = LinRange(0, 2π, Nᾱ)
#ᾱlist = LinRange(0, π/2, Nᾱ)
ψlist = LinRange(0.35, 0.65, Nκ)
θlist = LinRange(0, 2π, Nᾱ)
#%%
κ_res = zeros(Nκ, Nᾱ);
ᾱ_res = zeros(Nκ, Nᾱ);
α_res = zeros(Nκ, Nᾱ);
α2_res = zeros(Nκ, Nᾱ);
ψ_res = zeros(Nκ, Nᾱ);
θ_res = zeros(Nκ, Nᾱ);

for i in 1:Nκ, j in 1:Nᾱ
    ψ, θ, _ = MID.Mapping.isl_in_coords_to_tor(κlist[i], ᾱlist[j], ζval, isl)
    κ, ᾱ, _ = MID.Mapping.tor_coords_to_isl(ψlist[i], θlist[j], ζval, isl)
    ψ_res[i, j] = ψ
    θ_res[i, j] = θ
    κ_res[i, j] = κ
    ᾱ_res[i, j] = ᾱ
    α2_res[i, j] = mod((θlist[j] - ζval / isl.q0) * isl.m0, 2π)
    α_res[i, j] = 2/isl.m0 * asin(sqrt(κlist[i])*Elliptic.Jacobi.sn(2*Elliptic.K(κlist[i]) * ᾱlist[j]/π, κlist[i]))
end
#%%
#perhaps these are problemos now?
#looks much better for (1, -1) island, but probably hinting at a problemo!
#don't love the jagged edges we are still getting
#looks like this wasn't ever working, even for the radial case.
#unsure if that matters. But should fix so we can test both directions of mapping properly
contourf(ᾱlist, κlist, ψ_res)
contourf(ᾱlist, κlist, θ_res)
contourf(ᾱlist, κlist, α_res)
#%%
contourf(θlist, ψlist, κ_res)
contourf(θlist, ψlist, ᾱ_res)
contourf(θlist, ψlist, α2_res)

