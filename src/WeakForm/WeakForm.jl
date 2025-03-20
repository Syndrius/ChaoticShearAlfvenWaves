"""

Computes the weak form of our governing equation. This involves computing W and I for a given grid point.
"""
module WeakForm

using MID.Geometry
using MID.MagneticField
using MID.Structures
using MID.QFM

using LinearAlgebra

#struct storing tempory matrices used for memory efficient weak form
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

include("Tl.jl") #name subject to change
include("Tj.jl") #name subject to change
include("W.jl")
include("I.jl")



"""
    W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, met::MetT, B::BFieldT, prob::ProblemT, r:: Array{Float64}, θ, ζ::AbstractRange)

Computes the two matrices W and I based on the weak form of the SAW governing equation.
Solving generalised eigenvalue problem Wϕ = ω^2Iϕ
"""
function W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, met::MetT, B::BFieldT, prob::ProblemT, r::Array{Float64}, θ::AbstractArray, ζ::AbstractArray, tm::TM)

    
    #compute the density.
    n = prob.dens.(r) :: Array{Float64}
    #TODO
    ωcap2 = ω_cap2.(r) :: Array{Float64}

    for k=1:1:length(ζ), j=1:1:length(θ), i=1:1:length(r)

        #compute the metric
        prob.compute_met(met, r[i], θ[j], ζ[k], prob.geo.R0)

        #compute the magnetic field.
        #everywhere else has met then B...
        compute_B!(B, met, prob.q, prob.isl, prob.isl2, r[i], θ[j], ζ[k])

        #computes the matrix D.
        compute_D!(met, B, tm.D)

        #compute the W matrix
        @views compute_W!(W[:, :, i, j, k], met, B, n[i], ωcap2[i], tm)

        #compute the I matrix
        @views compute_I!(I[:, :, i, j, k], met, B, n[i], prob.flr, tm.D, tm.F)

    end

    

end



"""
    W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, met::MetT, B::BFieldT, prob::ProblemT, r:: Array{Float64}, θ, ζ::AbstractRange)

Computes the two matrices W and I based on the weak form of the SAW governing equation for the case with island coordinates.
In this case a specific metric and q-profile are used and the current nerm is not included.
"""
function W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, met::MetT, B::BFieldT, prob::IslProblemT, r::Array{Float64}, θ::AbstractArray, ζ::AbstractArray, tm::TM)
    
    #compute the density.
    n = prob.dens.(r) :: Array{Float64}
    #wcap is unused in this case, perhaps density as well.


    for k=1:1:length(ζ), j=1:1:length(θ), i=1:1:length(r)

        #compute the metric
        island_metric!(met, r[i], θ[j], ζ[k], prob.geo.R0, prob.isl)

        #compute the magnetic field.
        compute_B_isl!(B, met, prob.isl, r[i], θ[j], ζ[k])

        #computes the matrix D.
        compute_D!(met, B, tm.D)

        #compute the W matrix
        @views compute_isl_W!(W[:, :, i, j, k], met, B, tm)

        #compute the I matrix
        @views compute_I!(I[:, :, i, j, k], met, B, n[i], prob.flr, tm.D, tm.F)

    end

end



#is W and I a stupid name?
#this kind of needs to be here because of the qfm strutcs,
#ideally this would not be the case.
#I guess we could store them in the structure module???
#seems sub-optimal.
#but this is not a good place to keep this function
#special case for qfm surfaces
function W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, tor_met::MetT, tor_B::BFieldT, qfm_met::MetT, qfm_B::BFieldT,   prob::ProblemT, s::Array{Float64}, ϑ::AbstractArray, φ::AbstractArray, tm::TM, surfs::SurfaceITPT, CT::CoordTsfmT)

    #compute the density.
    n = prob.dens.(s) :: Array{Float64}
    #TODO
    #ωcap2 = ω_cap2.(r) :: Array{Float64}
    ωcap2 = zeros(length(s))

    for k=1:1:length(φ), j=1:1:length(ϑ), i=1:1:length(s)

        #so I am unsure if the new coords are ever actually used???
        coord_transform!(s[i], ϑ[j], φ[k], CT, surfs)
        #so the inv is fkn wopping.
        #display(CT.JM)
        #display(CT.JM_inv)

        #compute the original metric
        #using the computed values of (r, θ, ζ)
        prob.compute_met(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], prob.geo.R0)
        #and original B field.
        compute_B!(tor_B, tor_met, prob.q, prob.isl, prob.isl2, CT.coords[1], CT.coords[2], CT.coords[3]) 
        #transform the metric
        met_transform!(tor_met, qfm_met, CT)

        #display(qfm_met.gl)
        #transform the B field
        B_transform!(tor_B, qfm_B, qfm_met, CT)

        #display(tor_met.gl)
        #display(qfm_met.gl)


        #now we compute the weakform in the usual way.
        #computes the matrix D.
        WeakForm.compute_D!(qfm_met, qfm_B, tm.D)

        #compute the W matrix
        @views WeakForm.compute_W!(W[:, :, i, j, k], qfm_met, qfm_B, n[i], ωcap2[i], tm)

        #compute the I matrix
        @views WeakForm.compute_I!(I[:, :, i, j, k], qfm_met, qfm_B, n[i], prob.flr, tm.D, tm.F)

        #as a comparison!
        #does seem to be working in the normal case.
        #@views WeakForm.compute_W!(W[:, :, i, j, k], tor_met, tor_B, n[i], ωcap2[i], tm)

        #compute the I matrix
        #@views WeakForm.compute_I!(I[:, :, i, j, k], tor_met, tor_B, n[i], prob.flr, tm.D, tm.F)

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
    function compute_D!(met::MetT, B::BFieldT, D::Array{Float64, 2})

Function to compute D matrix, used by both W and I.
"""
function compute_D!(met::MetT, B::BFieldT, D::Array{Float64, 2})

    #D^ij = g^ij - b^i b^j
    @inbounds for j in 1:3
        @inbounds for i in 1:3

            D[i, j] = met.gu[i, j] - B.b[i]*B.b[j]
            
        end
    end

end


"""
    function init_tm()

Initialises structure storing the temporary matrices needed for the weak form.
"""
function init_tm()

    return TM(zeros(3, 9), zeros(3, 3), zeros(9, 3), zeros(9), zeros(3, 3), zeros(3, 3, 3), zeros(6))
end


end
