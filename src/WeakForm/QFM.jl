
"""
    weak_form!(P::Array{ComplexF64, 5}, Q::Array{ComplexF64, 5}, tor_B::BFieldT, tor_met::MetT, qfm_B::BFieldT, qfm_met::MetT, prob::ProblemT, s::Array{Float64}, ϑ::AbstractArray, ζ::AbstractArray, tm::TM, surfs::SurfaceITPT, CT::CoordTransformT, sd::TempSurfT)

Computes the two matrices P and Q based on the weak form of the SAW governing equation for the case with qfm surfaces.
The surfaces are used to convert the (s, ϑ, ζ) grid into (ψ, θ, φ) values, then the original metric and B are computed.
These are then transformed into the B and metric in (s, ϑ, ζ) coordinates so that the weakform is computed in terms of (s, ϑ, ζ).
"""
function weak_form!(P::Array{ComplexF64, 5}, Q::Array{ComplexF64, 5}, tor_B::BFieldT, tor_met::MetT, qfm_B::BFieldT, qfm_met::MetT, prob::ProblemT, s::Array{Float64}, ϑ::AbstractArray, ζ::AbstractArray, tm::TM, surfs::SurfaceITPT, CT::CoordTransformT, sd::TempSurfT)

    #compute the density.
    n = prob.dens.(s) :: Array{Float64}

    for k=1:1:length(ζ), j=1:1:length(ϑ), i=1:1:length(s)

        #compute the original coords, (ψ, θ, φ) and the jacobian matrix of the transformation.
        coord_transform!(s[i], ϑ[j], ζ[k], CT, surfs, sd)

        #compute the original metric
        #using the computed values of (ψ, θ, φ)
        prob.met(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], prob.geo.R0)

        #and original B field.
        compute_B!(tor_B, tor_met, prob.fields.q, prob.fields.isls, CT.coords[1], CT.coords[2], CT.coords[3]) 

        #transform the metric
        met_transform!(tor_met, qfm_met, CT)

        #transform the B field
        B_transform!(tor_B, qfm_B, qfm_met, CT)

        #now we compute the weakform in the usual way.

        #computes the matrix D.
        compute_D!(qfm_B, qfm_met, tm.D)

        #compute the W matrix
        @views compute_W!(W[:, :, i, j, k], qfm_B, qfm_met, n[i], tm)

        #compute the I matrix
        @views compute_I!(I[:, :, i, j, k], qfm_B, qfm_met, n[i], prob.flr, tm.D, tm.F)


    end

end

