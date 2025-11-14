#this is still v unclear
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
@kwdef struct RadialIslandT <: IslandT
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
    ψ0 :: Float64 = NaN
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
        isl = CoordIslandT(m0=m0, n0=n0, A=A, q0 = -m0/n0, qp=qp, ψ0 = ψ0, w=w)
    else
        isl = RadIslandT(m0=m0, n0=n0, A=A, q0 = -m0/n0, qp=qp, r0 = r0, w=w)
    end

    return isl

end


