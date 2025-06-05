#ideally this will one day have a qfm transform in here as well, although the inverse qfm transform may be cooked af.

function tor_to_qfm_f!(F::Array{Float64}, qfm_coords::Array{Float64}, tor_coords::Array{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

    coord_transform!(qfm_coords..., CT, surf_itp, sd)

    F .= CT.coords .- tor_coords
end

function tor_to_qfm_j!(J::Array{Float64}, qfm_coords::Array{Float64}, tor_coords::Array{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

    coord_transform!(qfm_coords..., CT, surf_itp, sd)

    J .= CT.JM
end

function tor_to_qfm_fj!(F::Array{Float64}, J::Array{Float64, 2}, qfm_coords::Array{Float64}, tor_coords::Array{Float64}, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

    coord_transform!(qfm_coords..., CT, surf_itp, sd)

    F .= CT.coords .- tor_coords
    J .= CT.JM
end

#there is no way this won't take a billion years
function tor_coords_to_qfm(r::Float64, θ::Float64, ζ::Float64, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)    

    f!(F, x) = tor_to_qfm_f!(F, x, [r, θ, ζ], CT, surf_itp, sd)
    j!(J, x) = tor_to_qfm_j!(J, x, [r, θ, ζ], CT, surf_itp, sd)
    fj!(F, J, x) = tor_to_qfm_fj!(F, J, x, [r, θ, ζ], CT, surf_itp, sd)

    #perhaps this should be defined earlier?
    x0 = [r, θ, ζ]
    F0 = similar(x0)


    #is 0.0 a terrible guess? seems to cause issues in the other case.
    df = OnceDifferentiable(f!, j!, fj!, x0, F0)

    sol = nlsolve(df, x0)

    s, ϑ, φ = sol.zero

    #think this is only needed for small grids but who knows
    #written in vector form to match sol.zero form.
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



#have decided that τ will be the new toroidal coordinates, even though it is unchanged.
#going to split the function up tbh. much easier.
function isl_in_coords_to_tor(κ::Float64, ᾱ::Float64, τ::Float64, isl::IslandT)

    #assumes κ < 1, throughout!
    sinβ = Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ)
    #α = mod(2/isl.m0 * asin(sqrt(κ)*sinβ), 2π)
    α = 2/isl.m0 * asin(sqrt(κ)*sinβ)

    α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ))

    θ = mod(α + τ/isl.q0, 2π)

    ζ = mod(τ, 2π)


    #abs herre prevents cases with ~-e-20
    r2diff = sqrt(abs(isl.w^2 * (κ - sin(isl.m0*α/2)^2)))

    #r2diff = sqrt(isl.w^2*(κ - sin(isl.m0*α /2)^2))

    #inside case, consider quadrants

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

#unsure where this belongs tbh
#computes the sepratrix in toroidal coodinates.
#note that we have a version of this in post-processing which just gives the min and max of the sepratrix.
#pretty stupid tbh
function compute_sepratrix(grids, isl)

    sep1 = zeros(grids.θ.N)
    sep2 = zeros(grids.θ.N)

    #idealy we can generalise this.
    #qp = isl.qp
    #q0 = isl.q0
    #r0 = isl.r0

    _, θgrid, _ = inst_grids(grids)

    #κ = 4/w^2(r^2/2 - r0^2/2)^2 + sin^@(m0 α / 2)
    # => κ = 1/w^2(r^2 - r0^2)^2 + sin^@(m0 α / 2)

    #this is going to assume ζ=0 for simplicity, should generalise though!
    for i in 1:grids.θ.N
        α = θgrid[i]
        #unsure why this needs to be isl.A + isl.A -> should be negative???
        #res = sqrt(abs(-2 * q0^2 / qp * (isl.A + isl.A * cos(isl.m0 * α))))

        #should define this in terms of width rather than all of this garbage.
        #res = sqrt(abs(16 * q0^2 * r0 * isl.A / qp * (1 -  sin(isl.m0 * α/2)^2)))

        #sep is defined for κ=1
        res = sqrt(isl.w^2 * (1 - sin(isl.m0*α/2)^2))

        #ψsep1[i] = sqrt(2*(-res + 0.125))
        #ψsep2[i] = sqrt(2*(res + 0.125))
        #ψsep1[i] = -sqrt(2*(res)) + 0.5# + 0.125))
        #ψsep2[i] = sqrt(2*(res)) + 0.5# + 0.125))
        sep1[i] = sqrt(-res + isl.r0^2)
        sep2[i] = sqrt(res + isl.r0^2)
    end

    return sep1, sep2
end
