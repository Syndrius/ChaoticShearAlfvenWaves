
"""
    separatrix(α::Float64, isl::RadialIslandT)

Computes the two radial values of the sepratrix for an input α.
"""
function separatrix(α::Float64, isl::RadialIslandT)

    r2diff = sqrt(isl.w^2*(1-sin(isl.m0*α/2)^2))

    return sqrt(-r2diff + isl.r0^2), sqrt(r2diff+isl.r0^2)
end

"""
    separatrix(α::Float64, isl::FluxIslandT)

Computes the two radial values of the sepratrix for an input α.
"""
function separatrix(α::Float64, isl::FluxIslandT)

    Δψ = sqrt(isl.w^2/4*(1-sin(isl.m0*α/2)^2))

    return -Δψ + isl.ψ0, Δψ + isl.ψ0
end


"""
    compute_sepratrix(θgrid::AbstractArray{Float64}, isl::IslandT, ζval::Float64=0.0)

Computes the radial values for of the sepratrix for each point in the θgrid.
"""
function compute_sepratrix(θgrid::AbstractArray{Float64}, isl::IslandT, ζval::Float64=0.0)

    sep1 = zeros(length(θgrid))
    sep2 = zeros(length(θgrid))


    for i in 1:length(θgrid)

        α = θgrid[i] + isl.n0/isl.m0 * ζval
        sep_min, sep_max = sepratrix(α, isl)

        sep1[i] = sep_min
        sep2[i] = sep_max
    end
    return sep1, sep2
end


"""
Fills in the extra values for the RadialIslandT, needed for the island q-profile.
"""
function inst_island(isl::RadialIslandT)

    q0 = -isl.m0/isl.n0

    qp = isl.qp

    r0 = isl.r0
    
    if isnan(isl.w)

        w = 4 * sqrt(q0^2*r0*isl.A / qp)
        A = isl.A
    else
        A = isl.w^2 / 16 * qp / (q0^2 * r0)
        w = isl.w
    end

    return RadialIslandT(isl.m0, isl.n0, A, q0, qp, r0, w)
end
    

"""
Fills in the extra values for the FluxIslandT, needed for the island q-profile.
"""
function inst_island(isl::FluxIslandT)

    q0 = -isl.m0/isl.n0

    qp = isl.qp

    ψ0 = isl.ψ0
    
    if isnan(isl.w)

        w = 4 * sqrt(q0^2*isl.A / qp)
        A = isl.A
    else
        A = isl.w^2 / 16 * qp / (q0^2)
        w = isl.w
    end

    return FluxIslandT(isl.m0, isl.n0, A, q0, qp, ψ0, w)
end
    
"""
Fills in the extra values for the CoordIslandT, needed for the island q-profile.
"""
function inst_island(isl::CoordIslandT)

    q0 = -isl.m0/isl.n0

    qp = isl.qp

    ψ0 = isl.ψ0
    
    if isnan(isl.w)

        w = 4 * sqrt(q0^2*isl.A / qp)
        A = isl.A
    else
        A = isl.w^2 / 16 * qp / (q0^2)
        w = isl.w
    end

    return CoordIslandT(isl.m0, isl.n0, A, isl.q0, qp, ψ0, isl.w)
end


"""
Converts an island defined for the geometric radius to toroidal flux.
"""
function convert_isl(isl::RadialIslandT)

    qp = isl.qp / isl.r0
    A = (isl.w/4)^2 * qp / isl.q0^2
    ψ0 = isl.r0^2 / 2
    
    return FluxIslandT(isl.m0, isl.n0, A, isl.q0, qp, ψ0, isl.w)
end


"""
Converts an island defined for the toroidal flux to geometric radius.
"""
function convert_isl(isl::FluxIslandT)

    r0 = sqrt(2 * isl.ψ0)

    qp = isl.qp * r0
    A = isl.w^2 / 16 * qp / (isl.q0^2 * r0)
    return RadialIslandT(isl.m0, isl.n0, A, isl.q0, qp, r0, isl.w)
end
    
