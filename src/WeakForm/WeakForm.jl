#TODO -> module needs to be cleaned up to match the other's and to not be quite so garbage.

module WeakForm

using MID.Geometry
using MID.MagneticField
using MID.Structures

using LinearAlgebra

export W_and_I!

include("W.jl")
include("I.jl")



"""
    W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, met::MetT, B::BFieldT, prob::ProblemT, r:: Array{Float64}, θ, ζ::AbstractRange)

Computes the two matrices W and I based on the weak form of the SAW governing equation.
Solving generalised eigenvalue problem Wϕ = ω^2Iϕ

### Args
- W::Array{ComplexF64, 5} - Matrix storing 9x9 values to be contracted with the derivtaives of ϕ, at each of the three coordinates.
- I::Array{ComplexF64, 5} - Matrix storing 9x9 values to be contracted with the derivtaives of ϕ, at each of the three coordinates.
- met::MetT - Struct containing the metric information.
- B::BfieldT - Struct containing the magnetic field information.
- prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
- r::Float64 - Radial coordinate, 0≤r≤1, minor radius is assumed 1.
- θ::Float64 - Poloidal angle, 0≤θ≤2π.
- ζ::Float64 - Toroidal angle, 0≤θ≤2π.
"""
function W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, met::MetT, B::BFieldT, prob::ProblemT, r::Array{Float64}, θ::AbstractArray, ζ::AbstractArray)

    #have removed the type for θ, needs to be a vector for zf, but a linrange for normal, need a better way of defining these types.
    
    #compute the density.
    n = prob.dens.(r) :: Array{Float64}
    #hard coded just for a test atm
    #dont actually need this!
    #dn = @. -5 * sech(8-10*r)^2

    for k=1:1:length(ζ), j=1:1:length(θ), i=1:1:length(r)

        #compute the metric
        prob.compute_met(met, r[i], θ[j], ζ[k], prob.geo.R0)

        #compute the magnetic field.
        compute_B!(B, met, prob.q, prob.isl, r[i], θ[j], ζ[k])
        #note that this is actually wrong, B[2, 2] was being overwritten with zero.
        #MagneticField.old_compute_B!(B, met, prob.q, prob.isl, r[i], θ[j], ζ[k])
        
        #comething about views here doesn't work, I assume it is do to with passing non-indexed things in ala met etc but who knows.
        #compute the W matrix
        #@views new_compute_W!(W[:, :, i, j, k], met, B)
        @views compute_W!(W[:, :, i, j, k], met, B)

        #views are giving us some warnings, may be better to just define the view of I/W not the entire @views thing.
        #compute the I matrix
        #@views new_compute_I!(I[:, :, i, j, k], met, B, n[i], prob.δ)
        @views compute_I!(I[:, :, i, j, k], met, B, n[i], prob.flr)
        #@views newest_compute_I!(I[:, :, i, j, k], met, B, n[i], prob.δ, dn[i], r[i])

    end

end

end