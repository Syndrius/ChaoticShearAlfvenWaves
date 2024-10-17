
#can we plot an 'island mode' in toroidal coordintaes.

using MID
using Plots#; plotlyjs()
using LaTeXStrings


function island_mode(κ, ᾱ, φ)

    #m=2 mode seems to have a `cleaner` fourier result.
    m = 2 
    n = 0

    if κ > 1 
        return 0
    end

    #approximate the quadratic looking global mode from Axel's paper.
    #-ve here dictates the orientation of peaks etc.
    #return (4*(κ - 0.5)^2 + 1)*exp(1im*m * ᾱ) #no toroidal component.

    #approximate a continuum mode.

    return exp(-(κ-0.7)^2/0.005) * exp(1im * m * ᾱ)

end


isl_grids = MID.Structures.MapGridsT(500, 2, 100, 10)

ϕ_isl = Array{ComplexF64}(undef, isl_grids.Nκ, isl_grids.Nᾱ, isl_grids.Nφ);

#κgrid, ᾱgrid, φgrid = inst_grids(isl_grids);

κgrid = LinRange(0, isl_grids.κmax, isl_grids.Nκ)
ᾱgrid = LinRange(0, 2π, isl_grids.Nᾱ)
φgrid = LinRange(0, 2π, isl_grids.Nφ)


for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, φ) in enumerate(φgrid)
    ϕ_isl[i, j, k] = island_mode(κ, ᾱ, φ)

end

#isl = IslandT(2, -1, 0.0006250000000000001, 2.0, 2.0, 0.5, 0.1)
isl = IslandT(2, -1, 5.625e-1, 2.0, 2.0, 0.5, 0.06)
tor_grids = MID.Structures.MapGridsT(800, 1, 400, 10)

#rgrid, θgrid, ζgrid = inst_grids(tor_grids)
rgrid = LinRange(0, tor_grids.κmax, tor_grids.Nκ)
θgrid = LinRange(0, 2π, tor_grids.Nᾱ)
ζgrid = LinRange(0, 2π, tor_grids.Nφ)

ϕ_tor, ϕfft_tor = isl_to_tor(tor_grids, ϕ_isl, isl_grids, isl);

heatmap(θgrid, rgrid, real.(ϕ_tor[:, :, 1]), levels=100, color=:turbo, ylimits=(0.4, 0.6), legend=false, xlabel=L"\theta", ylabel=L"r", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=25, xguidefontsize=25, xtickfontsize=10, ytickfontsize=10, dpi=600, xlimits=(0, 2π))


res = @. sqrt(isl.w^2 * (1 - sin(isl.m0 * θgrid / 2)^2));

seplow = @. sqrt(-res + isl.r0^2)
septop = @. sqrt(res + isl.r0^2)


scatter!(θgrid, seplow, color=:black, markersize=2, legend=false)
scatter!(θgrid, septop, color=:black, markersize=2, legend=false)

savefig("/scratch/y08/mt3516/aapps_data/isl_cont_mode.png")

#really highlights that the zero-zero mode is the most dominant...
potential_plot(ϕfft_tor, tor_grids)





using FFTW

ϕfft_tor_n0 =  Array{ComplexF64}(undef, tor_grids.Nκ, tor_grids.Nᾱ, tor_grids.Nφ);

ϕfft_tor_n0[:, :, :] = ϕfft_tor[:, :, :];


ϕfft_tor_n0[:, 1, 1] .= 0

potential_plot(ϕfft_tor_n0, tor_grids)

ϕ_tor_n0 = ifft(ϕfft_tor_n0, [2, 3]);

#probably not what was expected.
contourf(θgrid, rgrid, real.(ϕ_tor_n0[:, :, 1]), levels=100, color=:turbo)