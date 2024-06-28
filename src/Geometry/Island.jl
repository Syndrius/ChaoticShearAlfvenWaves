
"""
Struct storing the island parameters.

### Fields
- m0::Int64 - The poloidal mode number of the island chain.
- n0::Int64 - The toroidal mode number of the island chain.
- A::Float64 - The amplitude of the island. Requires more explaination...
"""
@kwdef struct IslandT
    m0 :: Int64 
    n0 :: Int64
    A :: Float64
end


"""
Struct storing the island parameters needed for computing the island continuum.

### Fields
- m0::Int64 - The poloidal mode number of the island chain.
- n0::Int64 - The toroidal mode number of the island chain.
- A::Float64 - The anplitude of the island #may need to relate this to width or explain it more
- q0::Float64 - value of q-profile at location of island.
- qp::Float64 - Derivative of q-profile at location of island.
- ψ0::Float64 - Radial location of island as a poloidal flux.
"""
@kwdef struct ContIslandT
    m0 :: Int64 
    n0 :: Int64
    A :: Float64
    q0 :: Float64
    qp :: Float64
    ψ0 :: Float64 #ideally change this to r if possible
end