
"""
    isl_in_coords_to_tor(κ::Float64, ᾱ::Float64, τ::Float64, isl::IslandT)

Maps island coordinates (κ, ᾱ, τ) into the equivalent toroidal coordinates (r, θ, φ), using the geometric radius.
Assumes that κ < 1.
"""
function isl_in_coords_to_tor(κ::Float64, ᾱ::Float64, τ::Float64, isl::RadialIslandT)

    #assumes κ < 1, throughout!
    sinβ = Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ)
    
    α = mod(2/isl.m0 * asin(sqrt(κ)*sinβ), 2π)
    θ = mod(α + τ/isl.q0, 2π)
    φ = mod(τ, 2π)

    #abs herre prevents cases with ~-e-20
    r2diff = sqrt(abs(isl.w^2 * (κ - sin(isl.m0*α/2)^2)))

    #split the solution into the different α quadrants
    if α < π/2
        #first quad, r>r0
        r = sqrt(+r2diff + isl.r0^2)
    elseif α < 3π/2
        #second or third quad r<r0
        r = sqrt(-r2diff+isl.r0^2)
    else 
        #fourth quad, r > r0
        r = sqrt(+r2diff+isl.r0^2)
    end
    return r, θ, φ
end


"""
    isl_in_coords_to_tor(κ::Float64, ᾱ::Float64, τ::Float64, isl::IslandT)

Maps island coordinates (κ, ᾱ, τ) into the equivalent toroidal coordinates (ψ, θ, φ), using the toroidal flux.
Assumes that κ < 1.
"""
function isl_in_coords_to_tor(κ::Float64, ᾱ::Float64, τ::Float64, isl::FluxIslandT)

    #assumes κ < 1, throughout!
    sinβ = Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ)

    α = mod(2/isl.m0 * asin(sqrt(κ)*sinβ), 2π/isl.m0)

    θ = mod(α + τ/isl.q0, 2π)
    φ = mod(τ, 2π)

    #abs herre prevents cases with ~-e-20
    Δψ = sqrt(abs(isl.w^2/4*(κ - sin(isl.m0*α/2)^2)))

    if ᾱ < π/2
        ψ = isl.ψ0 + Δψ 
    elseif ᾱ < 3π/2
        ψ = isl.ψ0 - Δψ 
    else 
        ψ = isl.ψ0 + Δψ 
    end
    return ψ, θ, φ
end

"""
    tor_coords_to_qfm(ψ::Float64, θ::Float64, ζ::Float64, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT) 

Maps toroidal coordinates (ψ, θ, φ) to qfm coordinates (s, ϑ, ζ). Requires root finding so can be slow.
"""
function tor_coords_to_qfm(ψ::Float64, θ::Float64, φ::Float64, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT) 

    f!(F, x) = tor_to_qfm_f!(F, x, [ψ, θ, φ], CT, surf_itp, sd)
    j!(J, x) = tor_to_qfm_j!(J, x, [ψ, θ, φ], CT, surf_itp, sd)
    fj!(F, J, x) = tor_to_qfm_fj!(F, J, x, [ψ, θ, φ], CT, surf_itp, sd)

    x0 = [ψ, θ, φ]
    F0 = similar(x0)

    #form needed by NLsolve.jl
    df = OnceDifferentiable(f!, j!, fj!, x0, F0)
    sol = nlsolve(df, x0)

    s, ϑ, ζ = sol.zero

    return (s, mod(ϑ, 2π), mod(ζ, 2π))
end


"""
    tor_coords_to_isl(r::Float64, θ::Float64, φ::Float64, isl::RadialIslandT)

Maps toroidal coordinates to island coordinates.
"""
function tor_coords_to_isl(r::Float64, θ::Float64, φ::Float64, isl::RadialIslandT)

    α = θ - φ / isl.q0

    κ = 1/isl.w^2 * (r^2 - isl.r0^2)^2 + sin(isl.m0 * α / 2)^2

    if κ < 1
        τ = asin(1/sqrt(κ) * abs(sin(isl.m0 * α / 2)))

        sinβ = 1/sqrt(κ) * sin(isl.m0 * α / 2)

        if r > isl.r0
            #i.e. top right qudrant.
            if mod(α * isl.m0, 2π) < π

                ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(τ, κ)

            else
                #top left quadrant, 4th by our clockwise notation.
                ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(-τ + 2π, κ)
            end
        else

            if mod(α * isl.m0, 2π) < π
                #bottom right quadrant.
                ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(-τ + π, κ)

            else
                #bottom left quadrant, 
                ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(τ + π, κ)
            end
        end
    else
        #this is still a work in progress.
        ᾱ = 0
        return 0, 0, 0
    end

    return κ, ᾱ, φ

end


"""
    tor_coords_to_isl(ψ::Float64, θ::Float64, φ::Float64, isl::IslandT)

Maps toroidal coordinates to island coordinates.
"""
function tor_coords_to_isl(ψ::Float64, θ::Float64, φ::Float64, isl::IslandT)

    α = θ - φ / isl.q0

    κ = 4/isl.w^2 * (ψ - isl.ψ0)^2 + sin(isl.m0 * α / 2)^2

    if κ < 1
        β = asin(1/sqrt(κ) * abs(sin(isl.m0 * α / 2)))

        if ψ == isl.ψ0
            return 0, 0, 0
        elseif ψ > isl.ψ0
            #i.e. top right qudrant.
            if mod(α * isl.m0, 2π) < π
                ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(β, κ)
            else
                #top left quadrant, 4th by our clockwise notation.
                ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(-β + 2π, κ)
            end
        else
            if mod(α * isl.m0, 2π) < π
                #bottom right quadrant.
                ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(-β + π, κ)
            else
                #bottom left quadrant, 
                ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(β + π, κ)
            end
        end
    else
        #this is still a work in progress.
        ᾱ = 0
        return 0, 0, 0
    end

    return κ, ᾱ, φ

end


"""
    tor_to_qfm_f!(F::Array{Float64}, qfm_coords::Array{Float64}, tor_coords::Array{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

Function used for root finding when mapping toroidal coords to qfm coords. Computes the residual.
"""
function tor_to_qfm_f!(F::Array{Float64}, qfm_coords::Array{Float64}, tor_coords::Array{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

    coord_transform!(qfm_coords..., CT, surf_itp, sd)

    F .= CT.coords .- tor_coords
end

"""
    tor_to_qfm_j!(J::Array{Float64, 2}, qfm_coords::Array{Float64}, tor_coords::Array{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

Function used for root finding when mapping toroidal coords to qfm coords. Computes the jacobian.
"""
function tor_to_qfm_j!(J::Array{Float64, 2}, qfm_coords::Array{Float64}, tor_coords::Array{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

    coord_transform!(qfm_coords..., CT, surf_itp, sd)

    J .= CT.JM
end

"""
    tor_to_qfm_fj!(F::Array{Float64}, J::Array{Float64, 2}, qfm_coords::Array{Float64}, tor_coords::Array{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

Function used for root finding when mapping toroidal coords to qfm coords. Computes the residual and the jacobian.
"""
function tor_to_qfm_fj!(F::Array{Float64}, J::Array{Float64, 2}, qfm_coords::Array{Float64}, tor_coords::Array{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

    coord_transform!(qfm_coords..., CT, surf_itp, sd)

    F .= CT.coords .- tor_coords
    J .= CT.JM
end


"""
    map_separatrix(ψmin::Array{Float64}, ψmax::Array{Float64}, θgrid::AbstractArray{Float64}, isl::IslandT, CT::CoordTransformT, surf_itp::QFM.SurfaceITPT, sd::TempSurfT)

Maps the island separatrix into qfm coordinates for deciding if solution is inside the island.
"""
function map_separatrix(ψmin::Array{Float64}, ψmax::Array{Float64}, θgrid::AbstractArray{Float64}, isl::IslandT, CT::CoordTransformT, surf_itp::QFM.SurfaceITPT, sd::TempSurfT)


    ψmin = zeros(length(θgrid))
    ψmax = zeros(length(θgrid))
    
    for (i, θ) in enumerate(θgrid)
        sep_min, sep_max = separatrix(θ, isl)
        ψmin[i] = sep_min
        ψmax[i] = sep_max
    end

    smin = zeros(length(ψmin))
    smax = zeros(length(ψmax))
    ϑsep = zeros(length(θgrid))

    for i in 1:length(θgrid)
        s, ϑ, _ = tor_coords_to_qfm(ψmin[i], θgrid[i], 0.0, CT, surf_itp, sd)
        smin[i] = s
        ϑsep[i] = ϑ 
        s, ϑ, _ = tor_coords_to_qfm(ψmax[i], θgrid[i], 0.0, CT, surf_itp, sd)
        smax[i] = s
    end
    return smin, smax, ϑsep

end
