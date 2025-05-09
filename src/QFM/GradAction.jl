"""
    action_grad_f!(δS::Array{Float64}, x::Array{Float64}, α::Float64, coefs::CoefficientsT, ipf::ActionGradFInputsT, ftd::FTDataT)

Computes the residuals of the extremized action equations.
"""
function action_grad_f!(δS::Array{Float64}, x::Array{Float64}, α::Float64, coefs::CoefficientsT, ipf::ActionGradFInputsT, ftd::FTDataT)
   
    unpack_coeffs!(x, coefs, ipf.aN+1)

    iota = ipf.b / ipf.a

    #computed area is just θ_0^c
    area = coefs.θcos[1]

    #ifft to get r, θ from the coefficients
    irfft1D!(ipf.r, coefs.rcos, coefs.rsin, ftd.irfft_1D_p, ftd.temp, ipf.aN, ipf.nfft)
    irfft1D!(ipf.θ, coefs.θcos, coefs.θsin, ftd.irfft_1D_p, ftd.temp, ipf.aN, ipf.nfft)

    #reflects that there is an additional term in the trial function for θ
    ipf.θ .+= iota .* ipf.ζ


    for i in 1:2*ipf.aN*ipf.nfft
        ipf.prob.met(ipf.met, ipf.r[i], ipf.θ[i], ipf.ζ[i], ipf.prob.geo.R0)

        compute_B!(ipf.B, ipf.met, ipf.prob.q, ipf.prob.isls, ipf.r[i], ipf.θ[i], ipf.ζ[i])

        #computes the `rhs' of the two equations.
        ipf.θdot.rhs[i] = ipf.B.B[2] / ipf.B.B[3]

        ipf.rdot.rhs[i] = ipf.B.B[1] / ipf.B.B[3] - coefs.ν[1] / (2*ipf.a*π*ipf.B.B[3] * ipf.met.J[1])
    end

    #fourier transform them for integration.
    rfft1D!(ipf.rdot.cn, ipf.rdot.sn, ipf.rdot.rhs, ftd.temp, ftd.rfft_1D_p, 2*ipf.nfft*ipf.aN)
    rfft1D!(ipf.θdot.cn, ipf.θdot.sn, ipf.θdot.rhs, ftd.temp, ftd.rfft_1D_p, 2*ipf.nfft*ipf.aN)


    #equations are structured so that first set comes from equation with rdot (i.e. δS/δθ) multipflied by cos
    #because cos(n=0) ≠ 0, we have aN + 1 equations.
    #second set comes from same equation multipflied by sin, giving aN equations
    #third set come from θdot eqaution (i.e. δS/δr) multipflied by cos
    #again, this has aN + 1 equations
    #fourth set is same equations with sin
    #giving aN equations
    #final set is a single equation for ν.
    #[rcos ; rsin[2:end] ; tcos ; tsin[2:end] ; [nv]]
    #δS has length 4*aN + 3.
    #[aN+1 ; aN ; aN+1 ; aN ; 1]

    #first aN+1 equations, from first equation integrated over cos
    δS[1:ipf.aN+1] = @. coefs.rsin * ipf.nlist / ipf.a - ipf.rdot.cn[1:ipf.aN+1]

    #second aN equations, from first equation integrated over sin
    #[2:end] as we dont need the rsin[1] coef, as it is zero
    δS[ipf.aN+2:2*ipf.aN+1] = @. (-coefs.rcos * ipf.nlist / ipf.a - ipf.rdot.sn[1:ipf.aN+1])[2:end]

    #third set of aN+1 equations, from first equation integrated over sin
    δS[2*ipf.aN+2:3*ipf.aN+2] = @. coefs.θsin * ipf.nlist / ipf.a - ipf.θdot.cn[1:ipf.aN+1]

    #(n=0) term for cos has additional iota term
    δS[2*ipf.aN+2] += iota

    #fourth set of aN equations, from first equation integrated over sin
    #[2:end] as we dont need the rsin[1] coef, as it is zero
    δS[3*ipf.aN+3:4*ipf.aN+2] = @. (-coefs.θcos * ipf.nlist / ipf.a - ipf.θdot.sn[1:ipf.aN+1])[2:end]

    #lagrange multipflier term
    δS[end] = area - α
end


"""
    action_grad_j!(JM::Array{Float64, 2}, x::Array{Float64}, α::Float64, coefs::CoefficientsT, ipf::ActionGradFInputsT, ftd::FTDataT, ipj::ActionGradJInputsT)

Computes the jacobian for the equations that extremize the action.
"""
function action_grad_j!(JM::Array{Float64, 2}, x::Array{Float64}, α::Float64, coefs::CoefficientsT, ipf::ActionGradFInputsT, ftd::FTDataT, ipj::ActionGradJInputsT)

    unpack_coeffs!(x, coefs, ipf.aN+1)

    iota = ipf.b / ipf.a

    #computed area is just θ_0^c
    area = coefs.θcos[1]

    #ifft to get r, θ from the coefficients
    irfft1D!(ipf.r, coefs.rcos, coefs.rsin, ftd.irfft_1D_p, ftd.temp, ipf.aN, ipf.nfft)
    irfft1D!(ipf.θ, coefs.θcos, coefs.θsin, ftd.irfft_1D_p, ftd.temp, ipf.aN, ipf.nfft)

    #reflects that there is an additional term in the trial function for θ
    ipf.θ .+= iota .* ipf.ζ


    for i in 1:2*ipf.aN*ipf.nfft

        ipf.prob.met(ipf.met, ipf.r[i], ipf.θ[i], ipf.ζ[i], ipf.prob.geo.R0)
        
        compute_B!(ipf.B, ipf.met, ipf.prob.q, ipf.prob.isls, ipf.r[i], ipf.θ[i], ipf.ζ[i])

        #derivatives of the `rhs' terms w.r.t r and θ
        rdot_rhs_dr = (ipf.B.dB[1, 1]/ipf.B.B[3] - ipf.B.B[1] * ipf.B.dB[3, 1] / ipf.B.B[3]^2 
                       + coefs.ν[1] * (ipf.B.dB[3, 1] / (2*ipf.a*π*ipf.B.B[3]^2*ipf.met.J[1]) 
                                        + ipf.met.dJ[1] / (2*ipf.a*π*ipf.B.B[3] * ipf.met.J[1]^2)))
        rdot_rhs_dθ = (ipf.B.dB[1, 2]/ipf.B.B[3] - ipf.B.B[1] * ipf.B.dB[3, 2] / ipf.B.B[3]^2 
                       + coefs.ν[1] * (ipf.B.dB[3, 2] / (2*ipf.a*π*ipf.B.B[3]^2*ipf.met.J[1]) 
                                        + ipf.met.dJ[2] / (2*ipf.a*π*ipf.B.B[3] * ipf.met.J[1]^2)))

        θdot_rhs_dr = ipf.B.dB[2, 1]/ipf.B.B[3] - ipf.B.B[2] * ipf.B.dB[3, 1] / ipf.B.B[3]^2

        θdot_rhs_dθ = ipf.B.dB[2, 2]/ipf.B.B[3] - ipf.B.B[2] * ipf.B.dB[3, 2] / ipf.B.B[3]^2

        #rhs of equation 1), derivative w.r.t ν
        ipj.ν_rhs[i] = 1 / (2*ipf.a*π * ipf.B.B[3] * ipf.met.J[1])
        
        #each derivative is multiplied by cos/sin(nζ/a)
        #as derivatives are taken with respect to the fourier coefficient
        #giving these terms as a chain rule
        for j in 1:ipf.aN+1
            #derivative w.r.t r^c
            ipj.drdot.drcn[i, j] = rdot_rhs_dr * ipj.cosnza[i, j]
            #derivative w.r.t r^s
            ipj.drdot.drsn[i, j] = rdot_rhs_dr * ipj.sinnza[i, j]
            #derivative w.r.t θ^c
            ipj.drdot.dθcn[i, j] = rdot_rhs_dθ * ipj.cosnza[i, j]
            #derivative w.r.t θ^s
            ipj.drdot.dθcn[i, j] = rdot_rhs_dθ * ipj.sinnza[i, j]

            #derivative w.r.t r^c
            ipj.dθdot.drcn[i, j] = θdot_rhs_dr * ipj.cosnza[i, j]
            #derivative w.r.t r^s
            ipj.dθdot.drsn[i, j] = θdot_rhs_dr * ipj.sinnza[i, j]
            #derivative w.r.t θ^c
            ipj.dθdot.dθcn[i, j] = θdot_rhs_dθ * ipj.cosnza[i, j]
            #derivative w.r.t θ^s
            ipj.dθdot.dθsn[i, j] = θdot_rhs_dθ * ipj.sinnza[i, j]
        end
    end

    #these are packed in the same way as our variables
    #this allows us to efficeintly fourier transform everything at once
    #then add to our Jacobian matrix all at once.
    ipj.drdot.drhs .= [ipj.drdot.drcn  ipj.drdot.drsn[:, 2:end]  ipj.drdot.dθcn  ipj.drdot.dθsn[:, 2:end]]
    ipj.dθdot.drhs .= [ipj.dθdot.drcn  ipj.dθdot.drsn[:, 2:end]  ipj.dθdot.dθcn  ipj.dθdot.dθsn[:, 2:end]]

    #take the fourier transform along the second dimension, to efficiently fourier trasnform everthing at once.
    rfft1D!(ipj.drdot.dcn, ipj.drdot.dsn, ipj.drdot.drhs, ftd.jm_temp, ftd.jmrfft_1D_p, 2*ipf.nfft*ipf.aN)
    rfft1D!(ipj.dθdot.dcn, ipj.dθdot.dsn, ipj.dθdot.drhs, ftd.jm_temp, ftd.jmrfft_1D_p, 2*ipf.nfft*ipf.aN)

    #does the ν rhs derivatives.
    #ideally, this would all be put into the grand drdot, but that doesn't seem to work.
    rfft1D!(ipj.ν_cn, ipj.ν_sn, ipj.ν_rhs, ftd.temp, ftd.rfft_1D_p, 2*ipf.nfft*ipf.aN)

    
    #note that JM is structured as JM[equation, derivative w.r.t var]
    #so JM[1, aN+2] is derivative of the first equation, w.r.t r^s_1
    JM .= 0.0

    #Derivtive of first set of equations, w.r.t every other coefficient.
    #end - 1 here as the new drdot doesn't include the ν terms.
    JM[1:ipf.aN+1, 1:end-1] .= -ipj.drdot.dcn[1:ipf.aN+1, 1:end]

    #we are excluding the n=0 contribution
    JM[CartesianIndex.(2:ipf.aN+1, ipf.aN+2:2*ipf.aN+1)] += ipf.nlist[2:end] / ipf.a

    #derivatives of second set of equations.
    JM[ipf.aN+2:2*ipf.aN+1, 1:end-1] .= -ipj.drdot.dsn[2:ipf.aN+1, 1:end]

    JM[CartesianIndex.(ipf.aN+2:2*ipf.aN+1, 2:ipf.aN+1)] += -ipf.nlist[2:end]/ipf.a

    #derivatives of 3rd set of equations
    JM[2*ipf.aN+2:3*ipf.aN+2, 1:end-1] .= -ipj.dθdot.dcn[1:ipf.aN+1, 1:end]

    JM[CartesianIndex.(2*ipf.aN+3:3*ipf.aN+2, 3*ipf.aN+3:4*ipf.aN+2)] += ipf.nlist[2:end]/ipf.a

    #derivatives of 4th set of equations
    JM[3*ipf.aN+3:4*ipf.aN+2, 1:end-1] .= -ipj.dθdot.dsn[2:ipf.aN+1, 1:end]

    JM[CartesianIndex.(3*ipf.aN+3:4*ipf.aN+2, 2*ipf.aN+3:3*ipf.aN+2)] += -ipf.nlist[2:end]/ipf.a

    #extra term from ν equation
    #this is d/dθ_0^c.
    JM[end, 2*ipf.aN+2] += 1.0

    #finally, derivatives w.r.t ν are added for the first and second set of equations.
    JM[1:ipf.aN+1, end] += ipj.ν_cn[1:ipf.aN+1]
    JM[ipf.aN+2:2*ipf.aN+1, end] += ipj.ν_sn[2:ipf.aN+1]
end



"""
    action_grad_fj!(δS::Array{Float64}, x::Array{Float64}, α::Float64, coefs::CoefficientsT, ipf::ActionGradFInputsT, ftd::FTDataT, ipj::ActionGradJInputsT)

Computes the residuals and jacobian of the extremized action equations, used for efficient solving by NLSolve.
"""
function action_grad_fj!(δS::Array{Float64}, JM::Array{Float64, 2}, x::Array{Float64}, α::Float64, coefs::CoefficientsT, ipf::ActionGradFInputsT, ftd::FTDataT, ipj::ActionGradJInputsT)
    unpack_coeffs!(x, coefs, ipf.aN+1)

    iota = ipf.b / ipf.a

    #computed area is just θ_0^c
    area = coefs.θcos[1]

    #ifft to get r, θ from the coefficients
    irfft1D!(ipf.r, coefs.rcos, coefs.rsin, ftd.irfft_1D_p, ftd.temp, ipf.aN, ipf.nfft)
    irfft1D!(ipf.θ, coefs.θcos, coefs.θsin, ftd.irfft_1D_p, ftd.temp, ipf.aN, ipf.nfft)

    #reflects that there is an additional term in the trial function for θ
    ipf.θ .+= iota .* ipf.ζ

    for i in 1:2*ipf.aN*ipf.nfft
        ipf.prob.met(ipf.met, ipf.r[i], ipf.θ[i], ipf.ζ[i], ipf.prob.geo.R0)

        #need B.dB now, so we just compute the full thing!
        compute_B!(ipf.B, ipf.met, ipf.prob.q, ipf.prob.isls, ipf.r[i], ipf.θ[i], ipf.ζ[i])

        #computes the `rhs' of the two equations.
        ipf.rdot.rhs[i] = ipf.B.B[1] / ipf.B.B[3] - coefs.ν[1] / (2*ipf.a*π*ipf.B.B[3] * ipf.met.J[1])
        ipf.θdot.rhs[i] = ipf.B.B[2] / ipf.B.B[3]

        #derivatives of the `rhs' terms w.r.t r and θ
        rdot_rhs_dr = (ipf.B.dB[1, 1]/ipf.B.B[3] - ipf.B.B[1] * ipf.B.dB[3, 1] / ipf.B.B[3]^2 
                       + coefs.ν[1] * (ipf.B.dB[3, 1] / (2*ipf.a*π*ipf.B.B[3]^2*ipf.met.J[1]) 
                                        + ipf.met.dJ[1] / (2*ipf.a*π*ipf.B.B[3] * ipf.met.J[1]^2)))
        rdot_rhs_dθ = (ipf.B.dB[1, 2]/ipf.B.B[3] - ipf.B.B[1] * ipf.B.dB[3, 2] / ipf.B.B[3]^2 
                       + coefs.ν[1] * (ipf.B.dB[3, 2] / (2*ipf.a*π*ipf.B.B[3]^2*ipf.met.J[1]) 
                                        + ipf.met.dJ[2] / (2*ipf.a*π*ipf.B.B[3] * ipf.met.J[1]^2)))

        θdot_rhs_dr = ipf.B.dB[2, 1]/ipf.B.B[3] - ipf.B.B[2] * ipf.B.dB[3, 1] / ipf.B.B[3]^2

        θdot_rhs_dθ = ipf.B.dB[2, 2]/ipf.B.B[3] - ipf.B.B[2] * ipf.B.dB[3, 2] / ipf.B.B[3]^2

        #rhs of equation 1), derivative w.r.t ν
        ipj.ν_rhs[i] = 1 / (2*ipf.a*π * ipf.B.B[3] * ipf.met.J[1])
        
        #each derivative is multiplied by cos/sin(nζ/a)
        #as derivatives are taken with respect to the fourier coefficient
        #giving these terms as a chain rule
        for j in 1:ipf.aN+1
            #derivative w.r.t r^c
            ipj.drdot.drcn[i, j] = rdot_rhs_dr * ipj.cosnza[i, j]
            #derivative w.r.t r^s
            ipj.drdot.drsn[i, j] = rdot_rhs_dr * ipj.sinnza[i, j]
            #derivative w.r.t θ^c
            ipj.drdot.dθcn[i, j] = rdot_rhs_dθ * ipj.cosnza[i, j]
            #derivative w.r.t θ^s
            ipj.drdot.dθcn[i, j] = rdot_rhs_dθ * ipj.sinnza[i, j]

            #derivative w.r.t r^c
            ipj.dθdot.drcn[i, j] = θdot_rhs_dr * ipj.cosnza[i, j]
            #derivative w.r.t r^s
            ipj.dθdot.drsn[i, j] = θdot_rhs_dr * ipj.sinnza[i, j]
            #derivative w.r.t θ^c
            ipj.dθdot.dθcn[i, j] = θdot_rhs_dθ * ipj.cosnza[i, j]
            #derivative w.r.t θ^s
            ipj.dθdot.dθsn[i, j] = θdot_rhs_dθ * ipj.sinnza[i, j]
        end
    end

    #fourier transform them for integration.
    rfft1D!(ipf.rdot.cn, ipf.rdot.sn, ipf.rdot.rhs, ftd.temp, ftd.rfft_1D_p, 2*ipf.nfft*ipf.aN)
    rfft1D!(ipf.θdot.cn, ipf.θdot.sn, ipf.θdot.rhs, ftd.temp, ftd.rfft_1D_p, 2*ipf.nfft*ipf.aN)


    #these are packed in the same way as our variables
    #this allows us to efficeintly fourier transform everything at once
    #then add to our Jacobian matrix all at once.
    ipj.drdot.drhs .= [ipj.drdot.drcn  ipj.drdot.drsn[:, 2:end]  ipj.drdot.dθcn  ipj.drdot.dθsn[:, 2:end]]
    ipj.dθdot.drhs .= [ipj.dθdot.drcn  ipj.dθdot.drsn[:, 2:end]  ipj.dθdot.dθcn  ipj.dθdot.dθsn[:, 2:end]]

    #take the fourier transform along the second dimension, to efficiently fourier trasnform everthing at once.
    rfft1D!(ipj.drdot.dcn, ipj.drdot.dsn, ipj.drdot.drhs, ftd.jm_temp, ftd.jmrfft_1D_p, 2*ipf.nfft*ipf.aN)
    rfft1D!(ipj.dθdot.dcn, ipj.dθdot.dsn, ipj.dθdot.drhs, ftd.jm_temp, ftd.jmrfft_1D_p, 2*ipf.nfft*ipf.aN)

    #does the ν rhs derivatives.
    #ideally, this would all be put into the grand drdot, but that doesn't seem to work.
    rfft1D!(ipj.ν_cn, ipj.ν_sn, ipj.ν_rhs, ftd.temp, ftd.rfft_1D_p, 2*ipf.nfft*ipf.aN)


    #Residuals
    #first aN+1 equations, from first equation integrated over cos
    δS[1:ipf.aN+1] = @. coefs.rsin * ipf.nlist / ipf.a - ipf.rdot.cn[1:ipf.aN+1]

    #second aN equations, from first equation integrated over sin
    #[2:end] as we dont need the rsin[1] coef, as it is zero
    δS[ipf.aN+2:2*ipf.aN+1] = @. (-coefs.rcos * ipf.nlist / ipf.a - ipf.rdot.sn[1:ipf.aN+1])[2:end]

    #third set of aN+1 equations, from first equation integrated over sin
    δS[2*ipf.aN+2:3*ipf.aN+2] = @. coefs.θsin * ipf.nlist / ipf.a - ipf.θdot.cn[1:ipf.aN+1]

    #(n=0) term for cos has additional iota term
    δS[2*ipf.aN+2] += iota

    #fourth set of aN equations, from first equation integrated over sin
    #[2:end] as we dont need the rsin[1] coef, as it is zero
    δS[3*ipf.aN+3:4*ipf.aN+2] = @. (-coefs.θcos * ipf.nlist / ipf.a - ipf.θdot.sn[1:ipf.aN+1])[2:end]

    #lagrange multipflier term
    δS[end] = area - α
    
    #Now Jacobian Matrix
    JM .= 0.0

    #Derivtive of first set of equations, w.r.t every other coefficient.
    #end - 1 here as the new drdot doesn't include the ν terms.
    JM[1:ipf.aN+1, 1:end-1] .= -ipj.drdot.dcn[1:ipf.aN+1, 1:end]

    #we are excluding the n=0 contribution
    JM[CartesianIndex.(2:ipf.aN+1, ipf.aN+2:2*ipf.aN+1)] += ipf.nlist[2:end] / ipf.a

    #derivatives of second set of equations.
    JM[ipf.aN+2:2*ipf.aN+1, 1:end-1] .= -ipj.drdot.dsn[2:ipf.aN+1, 1:end]

    JM[CartesianIndex.(ipf.aN+2:2*ipf.aN+1, 2:ipf.aN+1)] += -ipf.nlist[2:end]/ipf.a

    #derivatives of 3rd set of equations
    JM[2*ipf.aN+2:3*ipf.aN+2, 1:end-1] .= -ipj.dθdot.dcn[1:ipf.aN+1, 1:end]

    JM[CartesianIndex.(2*ipf.aN+3:3*ipf.aN+2, 3*ipf.aN+3:4*ipf.aN+2)] += ipf.nlist[2:end]/ipf.a

    #derivatives of 4th set of equations
    JM[3*ipf.aN+3:4*ipf.aN+2, 1:end-1] .= -ipj.dθdot.dsn[2:ipf.aN+1, 1:end]

    JM[CartesianIndex.(3*ipf.aN+3:4*ipf.aN+2, 2*ipf.aN+3:3*ipf.aN+2)] += -ipf.nlist[2:end]/ipf.a

    #extra term from ν equation
    #this is d/dθ_0^c.
    JM[end, 2*ipf.aN+2] += 1.0

    #finally, derivatives w.r.t ν are added for the first and second set of equations.
    JM[1:ipf.aN+1, end] += ipj.ν_cn[1:ipf.aN+1]
    JM[ipf.aN+2:2*ipf.aN+1, end] += ipj.ν_sn[2:ipf.aN+1]

end

