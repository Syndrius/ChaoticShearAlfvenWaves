

module Mapping

using ..Structures
using ..Geometry
using ..Indexing #needed for grid_id and basis_id, seems like they should be in basis, also makes this module a bit awkward
#we may actually want to move interpolation into post-processing, as it is theoretically possible that we want to interpolate in the future
#that will cook the dependencies though!
using ..QFM #probably going to cause problemos
using ..PostProcessing
using ..Io


using Elliptic
#using FFTW #we will do the fft in postprocessing! This module will just map the ϕ -> ϕ
#using BSplineKit #doesn't do 2d.
using Interpolations #pretty fkn stupid that we need multiple interpolations packages. #we can probably actually remove this now!
using JLD2
using Printf
using NLsolve #unfor that this is needed



#export map_isl_to_tor!
#export map_tor_to_isl!


include("CoordTransform.jl")
include("HermiteInterpolation.jl") 
include("Eigenfunctions.jl")
include("Spectrum.jl") 


#probably should be somewhere else
#need more of these as well.
#maybe there is a way to combine them, but I am not sure

#this should only be called once, so will be allocating
function qfm_to_tor_coord_map(rgrid::AbstractArray{Float64}, θgrid::AbstractArray{Float64}, ζgrid::AbstractArray{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

    coord_map = Array{Tuple{Float64, Float64, Float64}}(undef, length(rgrid), length(θgrid), length(ζgrid))

    for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)
        coord_map[i, j, k] = tor_coords_to_qfm(r, θ, ζ, CT, surf_itp, sd)
    end

    return coord_map
                                                        
end

function qfm_to_isl_coord_map(κgrid::AbstractArray{Float64}, ᾱgrid::AbstractArray{Float64}, τgrid::AbstractArray{Float64}, isl::IslandT, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

    coord_map = Array{Tuple{Float64, Float64, Float64}}(undef, length(κgrid), length(ᾱgrid), length(τgrid))

    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, τ) in enumerate(τgrid)
        r, θ, ζ = isl_in_coords_to_tor(κ, ᾱ, τ, isl)
        coord_map[i, j, k] = tor_coords_to_qfm(r, θ, ζ, CT, surf_itp, sd)
    end

    return coord_map
                                                        
end


function tor_to_isl_coord_map(κgrid::AbstractArray{Float64}, ᾱgrid::AbstractArray{Float64}, τgrid::AbstractArray{Float64}, isl::IslandT)

    coord_map = Array{Tuple{Float64, Float64, Float64}}(undef, length(κgrid), length(ᾱgrid), length(τgrid))

    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, τ) in enumerate(τgrid)
        coord_map[i, j, k] = isl_in_coords_to_tor(κ, ᾱ, τ, isl)
    end

    return coord_map
                                                        
end

end
