"""
Struct storing the eigenvalues and associated data.

### Fields
- ω::Array{ComplexF64} - Array storing the normalised eigenvalues.
- r::Array{Float64} - Radial location of the dominant mode for each eigenvalue.
- modelabs::Array{Tuple{Int64, Int64}} - (m, n) label of the dominant mode for each eigenvalue.
"""
@kwdef struct EvalsT 
    ω :: Array{ComplexF64}
    r :: Array{Float64}
    modelabs :: Array{Tuple{Int64, Int64}}
end


"""
    find_ind(evals::EvalsT, val)

Finds the index of an eigenvalue, useful for plotting.
"""
function find_ind(evals::EvalsT, val)

    return argmin(abs.(abs.(evals.ω) .- val))
end


"""
Generic version of above, should be somewhere else.
"""
function find_ind(grid::AbstractArray, val)

    return argmin(abs.(grid .- val))
end
