

module Mapping

using ..Structures
using ..Geometry
using ..Indexing #needed for grid_id and basis_id, seems like they should be in basis, also makes this module a bit awkward
#we may actually want to move interpolation into post-processing, as it is theoretically possible that we want to interpolate in the future
#that will cook the dependencies though!


using Elliptic
#using FFTW #we will do the fft in postprocessing! This module will just map the ϕ -> ϕ
#using BSplineKit #doesn't do 2d.
using Interpolations #pretty fkn stupid that we need multiple interpolations packages.
#using JLD2
#using Printf



export map_isl_to_tor!
export map_tor_to_isl!



include("CoordTransform.jl")
include("HermiteInterpolation.jl") 
#include("Spectrum.jl") 
#don't need the structures anymore, will just use a normal grid rather than a MapGrid.


#grids probbaly need to be fff for now, mainly because the fft extends the grids.
#island also needs to be instantiated!
#this function takes an eigenfunction computed using island coordinates, and maps it to toroidal coordinates.
function map_isl_to_tor(tor_grids::GridsT, ϕ_isl::Array{ComplexF64, 3}, isl_grids::GridsT, isl::IslandT)
    #TODO

    rgrid, θgrid, ζgrid = inst_grids(tor_grids)

    #may need to add the periodicity into this.
    κgrid, ᾱgrid, φgrid = inst_grids(isl_grids)

    #TODO
    #despite what Interpolations.jl says, we need to have the functions containing the periodic element
    #so fkn cooked
    #may just be worth getting the Hermite stuff working!
    isl_itp = interpolate((κgrid, ᾱgrid, φgrid), ϕ_isl, (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));

    #this will obvs be a disaster for v values far from the island!
    isl_ext = extrapolate(isl_itp, Periodic());

    ϕ_tor = Array{ComplexF64}(undef, length(rgrid), length(θgrid), length(ζgrid));
    for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)

        κ, ᾱ, φ = tor_coords_to_isl(r, θ, ζ, isl)

        #originally, was 0, unsure exactly why this was here, or why it is needed.
        if κ==10000 && ᾱ==0 && φ==0
            #for some reason this case is being mapped to 2?
            ϕ_tor[i, j, k] = 0
        else

            ϕ_tor[i, j, k] = isl_ext(κ, ᾱ, φ)
        end
    end

    return ϕ_tor#, fft(ϕ_tor, [2, 3])

end

#anoying af that the arrays are abstract. This is because they are step ranges or whatever
function map_tor_to_isl!(ϕ_isl::Array{ComplexF64, 3}, κgrid::Array{Float64}, ᾱgrid::AbstractArray{Float64}, τgrid::AbstractArray{Float64}, ϕ_tor::Array{ComplexF64, 4}, rgrid::Array{Float64}, θgrid::AbstractArray{Float64}, ζgrid::AbstractArray{Float64}, isl::IslandT)

    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, τ) in enumerate(τgrid)

        #find the equivalent coordinate in toroidal coordinates.
        r, θ, ζ = isl_in_coords_to_tor(κ, ᾱ, τ, isl)

        #now we interpolate the value of ϕ_tor at this location.
        ϕ_isl[i, j, k] = hermite_interpolation(r, θ, ζ, ϕ_tor, rgrid, θgrid, ζgrid)

    end
end


#this version uses interpolation for the ϕ_tor.
function map_tor_to_isl!(ϕ_isl::Array{ComplexF64, 3}, κgrid::Array{Float64}, ᾱgrid::AbstractArray{Float64}, τgrid::AbstractArray{Float64}, ϕ_tor::Array{ComplexF64, 3}, rgrid::Array{Float64}, θgrid::AbstractArray{Float64}, ζgrid::AbstractArray{Float64}, isl::IslandT)

    #the periodic part of this is cooked.
    #can be fixed by periodifying the grids and ϕ, i.e. adding the first element to the end.
    tor_itp = interpolate((rgrid, θgrid, ζgrid), ϕ_tor, (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));

    tor_ext = extrapolate(tor_itp, Periodic());

    #iterate through each point on our island grid
    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, τ) in enumerate(τgrid)

        #find the equivalent coordinate in toroidal coordinates.
        #fkn worst name so far.
        r, θ, ζ = isl_in_coords_to_tor(κ, ᾱ, τ, isl)

        #now we interpolate this coordinate using our functions.
        ϕ_isl[i, j, k] = tor_ext(r, θ, ζ)
        

    end
    
end

end
