
"""
    isl_in_coords_to_tor(κ::Float64, ᾱ::Float64, τ::Float64, isl::IslandT)

Maps island coordinates (κ, ᾱ, τ) into the equivalent toroidal coordinates (r, θ, ζ).
Assumes that κ < 1.
"""
function rad_isl_in_coords_to_tor(κ::Float64, ᾱ::Float64, τ::Float64, isl::IslandT)

    #assumes κ < 1, throughout!
    sinβ = Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ)
    #may need to make sure the domain of ᾱ is correct
    α = mod(2/isl.m0 * asin(sqrt(κ)*sinβ), 2π)
    #α = 2/isl.m0 * asin(sqrt(κ)*sinβ)

    #α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ))

    θ = mod(α + τ/isl.q0, 2π)

    ζ = mod(τ, 2π)


    #abs herre prevents cases with ~-e-20
    r2diff = sqrt(abs(isl.w^2 * (κ - sin(isl.m0*α/2)^2)))

    #r2diff = sqrt(isl.w^2*(κ - sin(isl.m0*α /2)^2))

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
    return r, θ, ζ
end


function isl_in_coords_to_tor(κ::Float64, ᾱ::Float64, τ::Float64, isl::IslandT)

    #assumes κ < 1, throughout!
    sinβ = Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ)
    #may need to make sure the domain of ᾱ is correct
    α = mod(2/isl.m0 * asin(sqrt(κ)*sinβ), 2π)
    #α = 2/isl.m0 * asin(sqrt(κ)*sinβ)

    #α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ))

    θ = mod(α + τ/isl.q0, 2π)

    ζ = mod(τ, 2π)


    #abs herre prevents cases with ~-e-20
    r2diff = sqrt(abs(isl.w^2 * (κ - sin(isl.m0*α/2)^2)))

    Δψ = sqrt(abs(isl.w^2/4 * (κ - sin(isl.m0*α/2)^2)))

    #r2diff = sqrt(isl.w^2*(κ - sin(isl.m0*α /2)^2))

    #split the solution into the different α quadrants
    if α < π/2
        #first quad, r>r0
        #r = sqrt(+r2diff + isl.r0^2)
        ψ = Δψ + isl.ψ0
    elseif α < 3π/2
        #second or third quad r<r0
        #r = sqrt(-r2diff+isl.r0^2)
        ψ = -Δψ + isl.ψ0
    else 
        #fourth quad, r > r0
        #r = sqrt(+r2diff+isl.r0^2)
        ψ = Δψ + isl.ψ0
    end
    return ψ, θ, ζ
end

"""
    tor_coords_to_qfm(r::Float64, θ::Float64, ζ::Float64, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT) 

Maps toroidal coordinates (r, θ, ζ) to qfm coordinates (s, ϑ, φ). Requires root finding so can be slow.
"""
function tor_coords_to_qfm(r::Float64, θ::Float64, ζ::Float64, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT) 

    f!(F, x) = tor_to_qfm_f!(F, x, [r, θ, ζ], CT, surf_itp, sd)
    j!(J, x) = tor_to_qfm_j!(J, x, [r, θ, ζ], CT, surf_itp, sd)
    fj!(F, J, x) = tor_to_qfm_fj!(F, J, x, [r, θ, ζ], CT, surf_itp, sd)

    #perhaps this should be defined earlier?
    x0 = [r, θ, ζ]
    F0 = similar(x0)

    #form needed by NLsolve.jl
    df = OnceDifferentiable(f!, j!, fj!, x0, F0)
    sol = nlsolve(df, x0)

    s, ϑ, φ = sol.zero

    return (s, mod(ϑ, 2π), mod(φ, 2π))
end


#TODO
#not as interested in this mapping, but ideally we would get this to work.
#the outside stuff is going to be v tricky!
function tor_coords_to_isl(r::Float64, θ::Float64, ζ::Float64, isl::IslandT)

    #unsure on the domain of this, may be 2πm_0? or 2π/m0
    α = θ - ζ / isl.q0

    κ = 1/isl.w^2 * (r^2 - isl.r0^2)^2 + sin(isl.m0 * α / 2)^2

    if κ < 1
        #so this is at least partially wrong...
        τ = asin(1/sqrt(κ) * abs(sin(isl.m0 * α / 2)))

        sinβ = 1/sqrt(κ) * sin(isl.m0 * α / 2)

        if r > isl.r0
            #i.e. top right qudrant.
            #exact expressions need more explaining.
            #taken from restricted_mapping.jl, which uses abs to certify the correct quadrant,
            #may need to do that here as the 1/2 in the sin will probbaly cook things...
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


        #ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(τ, κ)
    else
        #this is still a work in progress.
        ᾱ = 0
        #without this we get some very surprising results, in that they match our real results pretty fkn well.
        return 0, 0, 0
    end

    return κ, ᾱ, ζ

end




function isl_out_coords_to_tor(κ::Float64, ᾱ::Float64, τ::Float64, isl::IslandT)

    γ = Elliptic.am(isl.m0 * Elliptic.K(1/κ) * ᾱ / π, 1/κ)
    α = 2 * γ / isl.m0

    θ = mod(α + τ/isl.q0, 2π)
    ζ = mod(τ, 2π)

    #assumes κ >= 1, also returns positive and negative r.
    r2diff = isl.w * sqrt(κ - sin(γ)^2)

    rp = sqrt(+r2diff + isl.r0^2)
    rm = sqrt(-r2diff + isl.r0^2)

    return rp, rm, θ, ζ

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
    map_sepratrix(rmin::Array{Float64}, rmax::Array{Float64}, θgrid::AbstractArray{Float64}, isl::IslandT, CT::CoordTransformT, surf_itp::QFM.SurfaceITPT, sd::TempSurfT)

Maps the sepratrix from toroidal coordinates to qfm coordinates. Used to determine island modes.
"""
function map_sepratrix(rmin::Array{Float64}, rmax::Array{Float64}, θgrid::AbstractArray{Float64}, isl::IslandT, CT::CoordTransformT, surf_itp::QFM.SurfaceITPT, sd::TempSurfT)


    #TODO, Probably actually want to plot this for hek sake.
    #should all be the same length
    smin = zeros(length(rmin))
    smax = zeros(length(rmax))
    ϑsep = zeros(length(θgrid))

    for i in 1:length(θgrid)
        #this is probbaly a wildy inefficient way of doing this
        #but it is only done ones so who cares.
        s, ϑ, _ =  Mapping.tor_coords_to_qfm(rmin[i], θgrid[i], 0.0, CT, surf_itp, sd)
        smin[i] = s
        ϑsep[i] = ϑ #will this be the same in both cases? probably...
        s, ϑ, _ = Mapping.tor_coords_to_qfm(rmax[i], θgrid[i], 0.0, CT, surf_itp, sd)
        smax[i] = s
    end


    return smin, smax, ϑsep

end
