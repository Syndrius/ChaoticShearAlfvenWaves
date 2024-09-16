

#handles coordinate mapping between island and non-island coords.
#think this at least mildy works,
#but this is awful
#requries a rewrite of island struct I think
#and redo of how island stuff is integrated into this code.
module Mapping

using Elliptic
using FFTW

using MID.Structures
#this is needed for contislandT, need to fix island stuff!
using MID.Geometry

include("Interpolation.jl")
include("ToroidalToIsland.jl")


export tor_to_isl


#this maps a toroidal function into an island function.
function tor_to_isl(Nκ, Nᾱ, Nφ, ϕ, grids::FFFGridsT, isl::ContIslandT)

    #note that the ϕ to be mapped must be fff and contain the deriv parts.
    #may implement others some day.

    #first we define the equivalent island grids

    #not sure how we are going to deal with κ max. may just need to be a parameter.
    #2 is completly arbitrary lol.
    κgrid = LinRange(0, 2, Nκ);
    #κgrid2 = LinRange(0, 1.5, Nκ2);
    ᾱgrid = LinRange(0, 2π, Nᾱ+1)[1:end-1];
    φgrid = LinRange(0, 2π, Nφ+1)[1:end-1];

    ϕ_isl = zeros(ComplexF64, Nκ, Nᾱ, Nφ);

    #iterate through each point on our island grid
    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, φ) in enumerate(φgrid)

        #find the equivalent coordinate in toroidal coordinates.
        #fkn worst name so far.
        r, θ, ζ = coords_isl_to_tor(κ, ᾱ, φ, isl)

        #now we interpolate this coordinate using our functions.

        ϕ_isl[i, j, k] = hermite_interpolation(r, mod(θ, 2π), mod(ζ, 2π), ϕ, grids)

    end

    #perhaps we should fft this???
    
    return ϕ_isl, fft(ϕ_isl, [2, 3])
    

end

end