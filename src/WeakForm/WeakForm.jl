"""

Computes the weak form of our governing equation. This involves computing W and I for a given grid point.
"""
module WeakForm

using MID.Geometry
using MID.Equilibrium
using MID.Structures
using MID.QFM

using LinearAlgebra


"""
Struct storing Tempory Matrices used for memory efficient weak form.
"""
struct TM 
    C :: Array{Float64, 2}
    D :: Array{Float64, 2}
    T :: Array{Float64, 2}
    F :: Array{Float64}
    Γ :: Array{Float64, 2}
    dΓ :: Array{Float64, 3}
    K :: Array{Float64}
    function TM()
        new(zeros(3, 9), zeros(3, 3), zeros(9, 3), zeros(9), zeros(3, 3), zeros(3, 3, 3), zeros(6))
    end
end


export W_and_I!
export TM
export init_tm

include("Tl.jl") 
include("Tj.jl") 
include("W.jl")
include("I.jl")



"""
    W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, B::BFieldT, met::MetT, prob::ProblemT, r:: Array{Float64}, θ, ζ::AbstractRange)

Computes the two matrices W and I based on the weak form of the SAW governing equation.
Solving generalised eigenvalue problem Wϕ = ω^2Iϕ
"""
function W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, B::BFieldT, met::MetT, prob::ProblemT, r::Array{Float64}, θ::AbstractArray, ζ::AbstractArray, tm::TM)

    
    #compute the density.
    n = prob.dens.(r) :: Array{Float64}
    #TODO
    ωcap2 = ω_cap2.(r) :: Array{Float64}

    for k=1:1:length(ζ), j=1:1:length(θ), i=1:1:length(r)

        #compute the metric
        prob.met(met, r[i], θ[j], ζ[k], prob.geo.R0)

        #compute the magnetic field.
        compute_B!(B, met, prob.q, prob.isls, r[i], θ[j], ζ[k])

        #computes the matrix D.
        compute_D!(B, met, tm.D)

        #compute the W matrix
        @views compute_W!(W[:, :, i, j, k], B, met, n[i], ωcap2[i], tm)

        #compute the I matrix
        @views compute_I!(I[:, :, i, j, k], B, met, n[i], prob.flr, tm.D, tm.F)

    end

end



"""
    W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, B::BFieldT, met::MetT, prob::IslProblemT, κ::Array{Float64}, ᾱ::AbstractArray, ζ::AbstractArray, tm::TM)

Computes the two matrices W and I based on the weak form of the SAW governing equation for the case with island coordinates.
In this case a specific metric and q-profile are used and the current nerm is not included.
"""
function W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, B::BFieldT, met::MetT, prob::IslProblemT, κ::Array{Float64}, ᾱ::AbstractArray, ζ::AbstractArray, tm::TM)
    
    #TODO, no longer working with island as an array
    #compute the density.
    #density is probably assumed to be flat over the island
    #unsure if we want this to be an option
    n = prob.dens.(κ) :: Array{Float64}
    #wcap is unused in this case, perhaps density as well.


    for k=1:1:length(ζ), j=1:1:length(ᾱ), i=1:1:length(κ)

        #compute the metric
        island_metric!(met, κ[i], ᾱ[j], ζ[k], prob.geo.R0, prob.isl)

        #compute the magnetic field.
        compute_B_isl!(B, met, prob.isl, κ[i], ᾱ[j], ζ[k])

        #computes the matrix D.
        compute_D!(B, met, tm.D)

        #compute the W matrix
        @views compute_isl_W!(W[:, :, i, j, k], B, met, tm)

        #compute the I matrix
        @views compute_I!(I[:, :, i, j, k], B, met, n[i], prob.flr, tm.D, tm.F)

    end

end


"""
    W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, tor_B::BFieldT, tor_met::MetT, qfm_B::BFieldT, qfm_met::MetT, prob::ProblemT, s::Array{Float64}, ϑ::AbstractArray, φ::AbstractArray, tm::TM, surfs::SurfaceITPT, CT::CoordTsfmT, sd::TempSurfT)

Computes the two matrices W and I based on the weak form of the SAW governing equation for the case with qfm surfaces.
The surfaces are used to convert the (s, ϑ, φ) grid into (r, θ, ζ) values, then the original metric and B are computed.
These are then transformed into the B and metric in (s, ϑ, φ) coordinates so that the weakform is computed in terms of (s, ϑ, φ).
"""
function W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, tor_B::BFieldT, tor_met::MetT, qfm_B::BFieldT, qfm_met::MetT, prob::ProblemT, s::Array{Float64}, ϑ::AbstractArray, φ::AbstractArray, tm::TM, surfs::SurfaceITPT, CT::CoordTsfmT, sd::TempSurfT)

    #compute the density.
    n = prob.dens.(s) :: Array{Float64}
    #TODO
    #ωcap2 = ω_cap2.(r) :: Array{Float64}
    ωcap2 = zeros(length(s))

    for k=1:1:length(φ), j=1:1:length(ϑ), i=1:1:length(s)

        #compute the original coords, (r, θ, ζ) and the jacobian matrix of the transformation.
        coord_transform!(s[i], ϑ[j], φ[k], CT, surfs, sd)

        #compute the original metric
        #using the computed values of (r, θ, ζ)
        prob.met(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], prob.geo.R0)

        #and original B field.
        compute_B!(tor_B, tor_met, prob.q, prob.isls, CT.coords[1], CT.coords[2], CT.coords[3]) 

        #transform the metric
        met_transform!(tor_met, qfm_met, CT)

        #transform the B field
        B_transform!(tor_B, qfm_B, qfm_met, CT)

        #now we compute the weakform in the usual way.

        #computes the matrix D.
        compute_D!(qfm_B, qfm_met, tm.D)

        #compute the W matrix
        @views WeakForm.compute_W!(W[:, :, i, j, k], qfm_B, qfm_met, n[i], ωcap2[i], tm)

        #compute the I matrix
        @views WeakForm.compute_I!(I[:, :, i, j, k], qfm_B, qfm_met, n[i], prob.flr, tm.D, tm.F)


    end

end


#TODO work in progress.
function ω_cap2(r::Float64)

    #this seems to be having a much larger effect for fff than ffs, interesting...
    β = 0.0000000000
    #stab in the dark lol.
    return β * (1-r)

end


"""
    function compute_D!(B::BFieldT, met::MetT, D::Array{Float64, 2})

Function to compute D matrix, used by both W and I.
"""
function compute_D!(B::BFieldT, met::MetT, D::Array{Float64, 2})

    #D^ij = g^ij - b^i b^j
    @inbounds for j in 1:3
        @inbounds for i in 1:3

            D[i, j] = met.gu[i, j] - B.b[i]*B.b[j]
            
        end
    end

end

end
