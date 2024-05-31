
"""
Struct storing the island parameters.

# Fields
m0::Int64 - The poloidal mode number of the island chain.
n0::Int64 - The toroidal mode number of the island chain.
A::Float64 - The anplitude of the island #may need to relate this to width or explain it more
"""
@kwdef struct IslandT
    m0 :: Int64 
    n0 :: Int64
    A :: Float64
end