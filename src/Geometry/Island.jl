
#maybe better if this was the default island case with m0, n0 and A. and the others are only used if needed.
abstract type IslandT end

"""
Struct storing the island parameters. Only m0, n0 and one of A or w are required.
Island takes form A*sin(m0*θ + n0*ζ) so m0 and n0 should have different sign.

### Fields
- m0::Int64 - The poloidal mode number of the island chain.
- n0::Int64 - The toroidal mode number of the island chain.
- A::Float64=NaN - The anplitude of the island. 
- q0::Float64=NaN - value of q-profile at location of island.
- qp::Float64=NaN - Derivative of q-profile at location of island.
- r0::Float64=NaN - Radial location of island.
- w::Float64=NaN - Island width in units of r^2/2.
"""
@kwdef struct RadIslandT <: IslandT
    m0 :: Int64 
    n0 :: Int64
    A :: Float64 = NaN
    q0 :: Float64 = NaN
    qp :: Float64 = NaN
    r0 :: Float64 = NaN
    w :: Float64 = NaN
end


@kwdef struct FluxIslandT <: IslandT
    m0 :: Int64 
    n0 :: Int64
    A :: Float64 = NaN
    q0 :: Float64 = NaN
    qp :: Float64 = NaN
    ψ0 :: Float64 = NaN
    w :: Float64 = NaN
end


#this is perhaps the only one that actually needs all of this info
#this is also not a very good name!
#we could also uniquely define this by r0?
#but I think our island equiv q actually uses all of this info
@kwdef struct CoordIslandT <: IslandT
    m0 :: Int64 
    n0 :: Int64
    A :: Float64 = NaN
    q0 :: Float64 = NaN
    qp :: Float64 = NaN
    r0 :: Float64 = NaN
    w :: Float64 = NaN
end

"""
    init_island(; w::Float64=NaN, m0::Int64, n0::Int64, qp::Float64=NaN, r0::Float64=NaN, A::Float64=NaN)

Initialises the island structure. Many of the extra parameters are filled in once the problem is fully defined.
"""
function init_island(; w::Float64=NaN, m0::Int64, n0::Int64, qp::Float64=NaN, r0::Float64=NaN, A::Float64=NaN, flux::Bool=false, ψ0::Float64=NaN, coords::Bool=false)

    if isnan(w) && isnan(A)
        display("Please define the island width or amplitude")
        return 0
    end
    
    
    if flux
        isl = FluxIslandT(m0=m0, n0=n0, A=A, q0 = -m0/n0, qp=qp, ψ0 = ψ0, w=w)
    elseif coords 
        isl = CoordIslandT(m0=m0, n0=n0, A=A, q0 = -m0/n0, qp=qp, r0 = r0, w=w)
    else
        isl = RadIslandT(m0=m0, n0=n0, A=A, q0 = -m0/n0, qp=qp, r0 = r0, w=w)
    end

    return isl

end


"""
    sepratrix(α::Float64, isl::RadIslandT)

Computes the two radial values of the sepratrix for an input α.
"""
function sepratrix(α::Float64, isl::RadIslandT)

    r2diff = sqrt(isl.w^2*(1-sin(isl.m0*α/2)^2))

    return sqrt(-r2diff + isl.r0^2), sqrt(r2diff+isl.r0^2)
end

"""
    sepratrix(α::Float64, isl::FluxIslandT)

Computes the two radial values of the sepratrix for an input α.
"""
function sepratrix(α::Float64, isl::FluxIslandT)

    Δψ = sqrt(isl.w^2/4*(1-sin(isl.m0*α/2)^2))

    return -Δψ + isl.ψ0, Δψ + isl.ψ0
end


"""
    compute_sepratrix(θgrid::AbstractArray{Float64}, isl::IslandT, ζval::Float64=0.0)

Computes the radial values for of the sepratrix for each point in the θgrid.
"""
function compute_sepratrix(θgrid::AbstractArray{Float64}, isl::RadIslandT, ζval::Float64=0.0)

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
    inst_island(isl::IslandT, q::Functions)

Fills in the remaining values of the island struct based on the q-profile.
"""
function inst_island(isl::IslandT, q::Function)

    #accounts for cases where q0 is not set.
    q0 = -isl.m0/isl.n0

    #display(isl.q0)
    #display(zero_q(0.0, isl, q))

    #creates a temport q-profile for finding the root.
    tmpq(r) = zero_q(r, q0, q)

    #display(tmpq(0.0))
    #display(tmpq(1.0))
    
    r0 = find_zero(tmpq, (0, 1), Bisection() )

    _,  qp = q(r0)

    #compute either the width or the amplitude, depending on which has been specified.
    if isnan(isl.w)

        #w = 4*sqrt(isl.A * isl.q0^2/qp)
        #note that this width is defined in units of r^2.
        #so width in terms of r^2, what will this even mean??
        #we may want a function that can actually compute the width in r^2/2 terms
        #actually define in terms of r^2/2 i.e. flux, so that this matches other cases.
        #so this width is in terms of the flux surfaces not the radius.
        #not actually a very useful parameter then
        #I guess it is kind of the width in from (0, 0.5)??
        w = 4 * sqrt(q0^2*r0*isl.A / qp)
        A = isl.A
    else
        #A = (isl.w / 4)^2 * qp / isl.q0^2
        A = isl.w^2 / 16 * qp / (q0^2 * r0)
        w = isl.w
    end
    
    #not ideal to be creating a new island struct. Accessors.jl did not work for this.
    return IslandT(m0=isl.m0, n0=isl.n0, A=A, q0=q0, qp=qp, r0=r0, w=w)

end


"""
    zero_q(r, isl, q_prof)

Placeholder q-profile for root finding.
"""
function zero_q(r, q0, q_prof)

    q, _ = q_prof(r)

    #display(q)

    return q - q0
end


"""
    inst_island(isl::IslandT)

Case for island coords where most information must be predefined.
"""
function inst_island(isl::RadIslandT)

    q0 = -isl.m0/isl.n0

    qp = isl.qp

    r0 = isl.r0
    
    if isnan(isl.w)

        w = 4 * sqrt(q0^2*r0*isl.A / qp)
        A = isl.A
    else
        #A = (isl.w / 4)^2 * qp / isl.q0^2
        A = isl.w^2 / 16 * qp / (q0^2 * r0)
        w = isl.w
    end

    return RadIslandT(m0=isl.m0, n0=isl.n0, A=A, q0=q0, qp=qp, r0=r0, w=w)
end
    
    
function inst_island(isl::CoordIslandT)

    q0 = -isl.m0/isl.n0

    qp = isl.qp

    r0 = isl.r0
    
    if isnan(isl.w)

        w = 4 * sqrt(q0^2*r0*isl.A / qp)
        A = isl.A
    else
        #A = (isl.w / 4)^2 * qp / isl.q0^2
        A = isl.w^2 / 16 * qp / (q0^2 * r0)
        w = isl.w
    end

    return CoordIslandT(m0=isl.m0, n0=isl.n0, A=A, q0=q0, qp=qp, r0=r0, w=w)
end

#TODO -> placeholder for the no_island case!
function convert_isl(isl::RadIslandT)
    return FluxIslandT(m0=1.0, n0=1.0, A=0.0)
end

    
