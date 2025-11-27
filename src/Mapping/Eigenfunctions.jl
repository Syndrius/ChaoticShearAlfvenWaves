

"""
    efunc_map!(ϕ_map::Array{ComplexF64, 3}, N1::Int64, N2::Int64, N3::Int64, ϕ::Array{ComplexF64, 4}, x1g::AbstractArray{Float64}, x2g::AbstractArray{Float64}, x3g::AbstractArray{Float64}, coord_map::Array{Tuple{Float64, Float64, Float64}})

Maps the eigenfunction to a new coordinate system using the coordinate map.
"""
function efunc_map!(ϕ_map::Array{ComplexF64, 3}, N1::Int64, N2::Int64, N3::Int64, ϕ::Array{ComplexF64, 4}, x1g::AbstractArray{Float64}, x2g::AbstractArray{Float64}, x3g::AbstractArray{Float64}, coord_map::Array{Tuple{Float64, Float64, Float64}})

    for i in 1:N1, j in 1:N2, k in 1:N3
        ϕ_map[i, j, k] = interpolation(coord_map[i, j, k]..., ϕ, x1g, x2g, x3g)
    end
end



"""
    qfm_to_tor_coord_map(rgrid::AbstractArray{Float64}, θgrid::AbstractArray{Float64}, ζgrid::AbstractArray{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

Creates a Array that maps qfm coordinates to toroidal coordinates. Allows for efficient interpolation.
"""
function qfm_to_tor_coord_map(ψgrid::AbstractArray{Float64}, θgrid::AbstractArray{Float64}, φgrid::AbstractArray{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

    coord_map = Array{Tuple{Float64, Float64, Float64}}(undef, length(ψgrid), length(θgrid), length(φgrid))

    for (i, ψ) in enumerate(ψgrid), (j, θ) in enumerate(θgrid), (k, φ) in enumerate(φgrid)
        coord_map[i, j, k] = tor_coords_to_qfm(ψ, θ, φ, CT, surf_itp, sd)
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
        ψ, θ, φ = isl_in_coords_to_tor(κ, ᾱ, τ, isl)
        coord_map[i, j, k] = tor_coords_to_qfm(ψ, θ, φ, CT, surf_itp, sd)
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


"""
    isl_to_tor_coord_map(rgrid::AbstractArray{Float64}, θgrid::AbstractArray{Float64}, ζgrid::AbstractArray{Float64}, isl::IslandT)

Creates an array that maps island coordinates to toroidal coordinates. Allows for efficient interpolation.
"""
function isl_to_tor_coord_map(ψgrid::AbstractArray{Float64}, θgrid::AbstractArray{Float64}, φgrid::AbstractArray{Float64}, isl::IslandT)

    coord_map = Array{Tuple{Float64, Float64, Float64}}(undef, length(ψgrid), length(θgrid), length(φgrid))

    for (i, ψ) in enumerate(ψgrid), (j, θ) in enumerate(θgrid), (k, φ) in enumerate(φgrid)
        coord_map[i, j, k] = tor_coords_to_isl(ψ, θ, φ, isl)
    end

    return coord_map
                                                        
end

