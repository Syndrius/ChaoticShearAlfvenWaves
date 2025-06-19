

"""
    efunc_map!(ϕ_map::Array{ComplexF64, 3}, N1::Int64, N2::Int64, N3::Int64, ϕ::Array{ComplexF64, 4}, x1g::AbstractArray{Float64}, x2g::AbstractArray{Float64}, x3g::AbstractArray{Float64}, coord_map::Array{Tuple{Float64, Float64, Float64}})

Maps the eigenfunction to a new coordinate system using the coordinate map.
"""
function efunc_map!(ϕ_map::Array{ComplexF64, 3}, N1::Int64, N2::Int64, N3::Int64, ϕ::Array{ComplexF64, 4}, x1g::AbstractArray{Float64}, x2g::AbstractArray{Float64}, x3g::AbstractArray{Float64}, coord_map::Array{Tuple{Float64, Float64, Float64}})

    for i in 1:N1, j in 1:N2, k in 1:N3
        ϕ_map[i, j, k] = hermite_interpolation(coord_map[i, j, k]..., ϕ, x1g, x2g, x3g)
    end
end



"""
    qfm_to_tor_coord_map(rgrid::AbstractArray{Float64}, θgrid::AbstractArray{Float64}, ζgrid::AbstractArray{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

Creates a Array that maps qfm coordinates to toroidal coordinates. Allows for efficient interpolation.
"""
function qfm_to_tor_coord_map(rgrid::AbstractArray{Float64}, θgrid::AbstractArray{Float64}, ζgrid::AbstractArray{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

    #note that there is nothing in here to distinguish radius from flux, just comes down to surfaces made I guess

    coord_map = Array{Tuple{Float64, Float64, Float64}}(undef, length(rgrid), length(θgrid), length(ζgrid))

    for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)
        coord_map[i, j, k] = tor_coords_to_qfm(r, θ, ζ, CT, surf_itp, sd)
    end

    return coord_map
                                                        
end

"""
    qfm_to_isl_coord_map(κgrid::AbstractArray{Float64}, ᾱgrid::AbstractArray{Float64}, τgrid::AbstractArray{Float64}, isl::IslandT, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

Creates an array that maps qfm coordinates to island coordinates. Allows for efficient interpolation.
"""
function qfm_to_isl_coord_map(κgrid::AbstractArray{Float64}, ᾱgrid::AbstractArray{Float64}, τgrid::AbstractArray{Float64}, isl::IslandT, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

    coord_map = Array{Tuple{Float64, Float64, Float64}}(undef, length(κgrid), length(ᾱgrid), length(τgrid))

    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, τ) in enumerate(τgrid)
        r, θ, ζ = isl_in_coords_to_tor(κ, ᾱ, τ, isl)
        coord_map[i, j, k] = tor_coords_to_qfm(r, θ, ζ, CT, surf_itp, sd)
    end

    return coord_map
                                                        
end


"""
    tor_to_isl_coord_map(κgrid::AbstractArray{Float64}, ᾱgrid::AbstractArray{Float64}, τgrid::AbstractArray{Float64}, isl::IslandT)

Creates an array that maps toroidal coordinates to island coordinates. Allows for efficient interpolation.
"""
function tor_to_isl_coord_map(κgrid::AbstractArray{Float64}, ᾱgrid::AbstractArray{Float64}, τgrid::AbstractArray{Float64}, isl::IslandT)

    coord_map = Array{Tuple{Float64, Float64, Float64}}(undef, length(κgrid), length(ᾱgrid), length(τgrid))

    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, τ) in enumerate(τgrid)
        coord_map[i, j, k] = isl_in_coords_to_tor(κ, ᾱ, τ, isl)
    end

    return coord_map
                                                        
end


function isl_to_tor_coord_map(rgrid::AbstractArray{Float64}, θgrid::AbstractArray{Float64}, ζgrid::AbstractArray{Float64}, isl::IslandT)

    coord_map = Array{Tuple{Float64, Float64, Float64}}(undef, length(rgrid), length(θgrid), length(ζgrid))

    for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)
        coord_map[i, j, k] = tor_coords_to_isl(r, θ, ζ, isl)
    end

    return coord_map
                                                        
end

