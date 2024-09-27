
"""
    island_W_and_I!(W::Array{ComplexF64, 2}, I::Array{ComplexF64, 2}, χ::Float64, θ̄grid::AbstractRange, ζgrid::AbstractRange, met::MetT, B::BFieldT, isl::IslandT, geo::GeoParamsT, sign::Int64)
    
Function to compute the two matrices, W and I, that make up the weak form in the island ocntinuum case,
essentially equation (23) of Qu and Hole, before Fourier transforming or acting the parallel gradient.

# Args
W::Array{ComplexF64, 2} - W matrix.
I::Array{ComplexF64, 2} - I matrix.
χ::Float64 - Energy value.
θ̄grid::Array{Float64} - Array of poloidal points.
ζgrid::Array{Float64} - Array of toroidal points.
met::MetT - Metric struct.
B::BFieldT - Magnetic field struct.
isl::IslandT - Struct storing island parameters.
geo::GeoParamsT - struct storing the geometrical parameters, i.e major radius.
sign::Int64 - Sign of the particles, ±1 for passing on either side, or 0 for trapped particles.
"""

function island_W_and_I!(W::Array{ComplexF64, 2}, I::Array{ComplexF64, 2}, χ::Float64, θ̄grid::AbstractRange, ζgrid::AbstractRange, met::MetT, B::BFieldT, isl::IslandT, geo::GeoParamsT, sign::Int64)

    for (i, θ̄) in enumerate(θ̄grid), (j, ζ) in enumerate(ζgrid)

        #first we determine what ᾱ is.
        if sign == 0
            ᾱ = θ̄
        else
            ᾱ = θ̄ - ζ/isl.q0
        end
        #now we determine our original coordinate values, for computing the metric and B.
        α = compute_α(χ, ᾱ, isl, sign)
        
        ψ = compute_ψ(χ, α, ᾱ, isl, sign)

        
        #theta is actually different for passing vs trapped,
        #but zetagrid differences handle that.
        θ = α + ζ/isl.q0

        flux_toroidal_metric!(met, ψ, θ, ζ, geo.R0)

        compute_island_B!(B, met, isl, ψ, α)

        ∇χ2 = compute_∇χ2(ψ, α, met, isl)


        W[i, j] = ∇χ2 / (met.J * B.mag_B^2)

        I[i, j] = ∇χ2 * met.J / B.mag_B^2
    end
        

end