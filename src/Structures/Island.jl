#abstract type for generic island
abstract type IslandT end

"""
Struct storing the island parameters. Only m0, n0 and one of A or w are required.
Island takes form A*sin(m0*θ + n0*ζ) so m0 and n0 should have different sign.
Different types are used for multiple dispatch based on radial coordinate.

### Fields
- m0::Int64 - The poloidal mode number of the island chain.
- n0::Int64 - The toroidal mode number of the island chain.
- A::Float64=NaN - The anplitude of the island. 
- q0::Float64=NaN - value of q-profile at location of island.
- qp::Float64=NaN - Derivative of q-profile at location of island.
- r0::Float64=NaN - Radial location of island.
- w::Float64=NaN - Island width in units of r^2/2.
"""

#geometric radius as radial coordinate
struct RadialIslandT <: IslandT
    m0 :: Int64 
    n0 :: Int64
    A :: Float64
    q0 :: Float64
    qp :: Float64
    r0 :: Float64
    w :: Float64
end

#toroidal flux as the radial coordinate
struct FluxIslandT <: IslandT
    m0 :: Int64 
    n0 :: Int64
    A :: Float64
    q0 :: Float64
    qp :: Float64
    ψ0 :: Float64
    w :: Float64
end


#Island for case where island straight field line coordinates are used.
struct CoordIslandT <: IslandT
    m0 :: Int64 
    n0 :: Int64
    A :: Float64
    q0 :: Float64
    qp :: Float64
    ψ0 :: Float64
    w :: Float64
end

const no_isl = FluxIslandT(0, 0, 0.0, 0.0, 0.0, 0.0, 0.0)

"""
    init_island(; w::Float64=NaN, m0::Int64, n0::Int64, qp::Float64=NaN, r0::Float64=NaN, A::Float64=NaN)

Initialises the island structure. Many of the extra parameters are filled in once the problem is fully defined.
"""
function init_island(type::Symbol=:ψ; w::Float64=NaN, m0::Int64, n0::Int64, qp::Float64=NaN, r0::Float64=NaN, A::Float64=NaN, ψ0::Float64=NaN)

    if isnan(w) && isnan(A)
        display("Please define the island width or amplitude")
        return 0
    end

    if type in [:r, :radial]
        isl = RadIslandT(m0, n0, A, -m0/n0, qp, r0, w)
    elseif type in [:κ, :island, :coords] 
        isl = CoordIslandT(m0, n0, A, -m0/n0, qp, ψ0, w)
    else
        isl = FluxIslandT(m0, n0, A, -m0/n0, qp, ψ0, w)
    end

    return isl

end

