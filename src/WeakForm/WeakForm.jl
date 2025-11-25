"""

Computes the weak form of our governing equation. This involves computing P and Q for a given grid point.
"""
module WeakForm

using ..Geometry 
using ..Fields
using ..Structures
using ..QFM


using LinearAlgebra


"""
Struct storing Tempory Matrices used for memory efficient weak form.
"""
struct TM 
    C :: Array{Float64, 2}
    D :: Array{Float64, 2}
    F :: Array{Float64}
    Γ :: Array{Float64, 2}
    dΓ :: Array{Float64, 3}
    K :: Array{Float64}
    function TM()
        new(zeros(3, 9), zeros(3, 3), zeros(9), zeros(3, 3), zeros(3, 3, 3), zeros(6))
    end
end


export weak_form!
export TM

#good
include("Tl.jl") 
#good
include("Tj.jl") 
#good
include("P.jl")
#good
include("Q.jl")
#good
include("QFM.jl")
#good
include("Problem.jl")

export init_problem, inst_problem


"""
    weak_form!(P::Array{ComplexF64, 5}, Q::Array{ComplexF64, 5}, B::BFieldT, met::MetT, prob::ProblemT, r:: Array{Float64}, θ, ζ::AbstractRange)

Computes the two matrices P and Q based on the weak form of the SAW governing equation.
Solving generalised eigenvalue problem PΦ = ω^2QΦ.
"""
function weak_form!(P::Array{ComplexF64, 5}, Q::Array{ComplexF64, 5}, B::BFieldT, met::MetT, prob::ProblemT, r::Array{Float64}, θ::AbstractArray, ζ::AbstractArray, tm::TM)

    
    #compute the density.
    n = prob.fields.dens.(r) :: Array{Float64}

    for k=1:1:length(ζ), j=1:1:length(θ), i=1:1:length(r)

        #compute the metric
        prob.geo.met(met, r[i], θ[j], ζ[k], prob.geo.R0)
        #display("or here")

        #compute the magnetic field.
        compute_B!(B, met, prob.fields.q, prob.fields.isls, r[i], θ[j], ζ[k])

        #computes the matrix D.
        compute_D!(B, met, tm.D)

        #compute the P matrix
        @views compute_P!(P[:, :, i, j, k], B, met, n[i], tm)

        #compute the Q matrix
        @views compute_Q!(Q[:, :, i, j, k], B, met, n[i], prob.flr, tm.D, tm.F)

    end

end


"""
    function compute_D!(B::BFieldT, met::MetT, D::Array{Float64, 2})

Function to compute D matrix, used by both P and Q.
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
