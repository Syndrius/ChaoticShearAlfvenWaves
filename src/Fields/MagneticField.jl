
"""
    compute_B!(B::BFieldT, met::MetT, q_prof::FunctionWrapper{Tuple{Float64, Float64}}, isls::Array{FluxIslandT}, ψ::Float64, x2::Float64, x3::Float64)

Function the fills the BFieldT struct based on an q-profile and islands. 
This function is different for different radial coordinates, dispatch is based on the type of island.
"""
function compute_B!(B::BFieldT, met::MetT, q_prof::FunctionWrapper{Tuple{Float64, Float64}}, isls::Array{FluxIslandT}, x1::Float64, x2::Float64, x3::Float64)

    q, dq = q_prof(x1)


    #loop through the islands
    B.B[1] = 0.0
    B.dB[1, :] .= 0.0
    for isl in isls
        B.B[1] += isl.A * isl.m0 * sin(isl.m0 * x2 + isl.n0 * x3)
        B.dB[1, 2] += isl.A * isl.m0^2 * cos(isl.m0 * x2 + isl.n0 * x3)
        B.dB[1, 3] += isl.A * isl.m0 * isl.n0 * cos(isl.m0 * x2 + isl.n0 * x3)
    end

    B.B[1] /= met.J[1]

    B.dB[1, 1] -= B.B[1] * met.dJ[1] / met.J[1]

    B.dB[1, 2] /= met.J[1]
    B.dB[1, 2] -= B.B[1] * met.dJ[2] / met.J[1]

    B.dB[1, 3] /= met.J[1]
    B.dB[1, 3] -= B.B[1] * met.dJ[3] / met.J[1] #typically dJ/dx3 = 0.
                


    B.B[2] = 1 / (met.J[1] * q) 
    B.B[3] = 1 / (met.J[1])

    B.dB[2, 1] = ( - (1 / q) * met.dJ[1] / met.J[1]^2
                    - 1 * dq /(met.J[1] * q^2) )
                    
    B.dB[2, 2] = - (1 / q ) * met.dJ[2] / met.J[1]^2

    B.dB[3, 2] = - (1 / q ) * met.dJ[3] / met.J[1]^2


    B.dB[3, 1] = ( - met.dJ[1]) / met.J[1]^2
    B.dB[3, 2] = - met.dJ[2] / met.J[1]^2
    B.dB[3, 3] = - met.dJ[3] / met.J[1]^2

    #note this also does B.db
    magnitude_B!(B, met)
    
    B.b[1] = B.B[1]/B.mag_B[1]
    B.b[2] = B.B[2]/B.mag_B[1]
    B.b[3] = B.B[3]/B.mag_B[1]
end


"""
This version is for the geometric radius.
"""
function compute_B!(B::BFieldT, met::MetT, q_prof::FunctionWrapper{Tuple{Float64, Float64}, Tuple{Float64}}, isls::Array{RadialIslandT}, x1::Float64, x2::Float64, x3::Float64)

    q, dq = q_prof(x1)


    B.B[1] = 0.0
    B.dB[1, :] .= 0.0
    for isl in isls
        B.B[1] += isl.A * isl.m0 * sin(isl.m0 * x2 + isl.n0 * x3)
        B.dB[1, 2] += isl.A * isl.m0^2 * cos(isl.m0 * x2 + isl.n0 * x3)
        B.dB[1, 3] += isl.A * isl.m0 * isl.n0 * cos(isl.m0 * x2 + isl.n0 * x3)
    end

    B.B[1] /= met.J[1]

    B.dB[1, 1] -= B.B[1] * met.dJ[1] / met.J[1]

    B.dB[1, 2] /= met.J[1]
    B.dB[1, 2] -= B.B[1] * met.dJ[2] / met.J[1]

    B.dB[1, 3] /= met.J[1]
    B.dB[1, 3] -= B.B[1] * met.dJ[3] / met.J[1] #typically dJ/dx3 = 0.

                
    B.B[2] = x1 / (met.J[1] * q) 
    B.B[3] = x1 / (met.J[1])


    B.dB[2, 1] = (1 / (met.J[1] * q) 
                    - (x1 / q) * met.dJ[1] / met.J[1]^2
                    - x1 * dq /(met.J[1] * q^2) )
                    
    B.dB[2, 2] = - (x1 / q ) * met.dJ[2] / met.J[1]^2
    

    B.dB[2, 3] = - (x1 / q ) * met.dJ[3] / met.J[1]^2 #typically dJ/dx3 = 0.

    B.dB[3, 1] = (met.J[1] - x1 * met.dJ[1]) / met.J[1]^2
    B.dB[3, 2] = -x1*met.dJ[2]/met.J[1]^2
    B.dB[3, 3] = -x1*met.dJ[3]/met.J[1]^2 #typically dJ/dx3 = 0.

    #note this also does B.db
    magnitude_B!(B, met)
    
    B.b[1] = B.B[1]/B.mag_B[1]
    B.b[2] = B.B[2]/B.mag_B[1]
    B.b[3] = B.B[3]/B.mag_B[1]
end


"""
    magnitude_B!(B::BFieldT, met::MetT)

Function that computes the magnitude of B and its derivative once the other fields of BFieldT have been computed.
"""
function magnitude_B!(B::BFieldT, met::MetT)

    B.mag_B .= 0.0
    B.dmag_B .= 0.0

    for i in 1:3, j in 1:3

        #note that this is actual |B|^2
        B.mag_B[1] += B.B[i] * B.B[j] * met.gl[i, j]

        for k in 1:3
            B.dmag_B[k] += (B.B[i] * B.B[j] * met.dgl[i, j, k] +
                        B.B[i] * met.gl[i, j] * B.dB[j, k] +
                        B.B[j] * met.gl[i, j] * B.dB[i, k])
        end
    end

    B.mag_B[1] = sqrt(B.mag_B[1])
    B.dmag_B[:] = @. 1/(2*B.mag_B[1]) * B.dmag_B

    for i in 1:3, j in 1:3
        B.db[j, i] = B.dB[j, i]/B.mag_B[1] - B.B[j]*B.dmag_B[i]/B.mag_B[1]^2
    end
end



"""
This version is for island coordinates.
"""
function compute_B!(B::BFieldT, met::MetT, q_prof::FunctionWrapper{Tuple{Float64, Float64}}, isls::Array{CoordIslandT}, κ::Float64, ᾱ::Float64, τ::Float64)


    isl = isls[1] #only a single island is cosnidered for this case.

    #should probably be q̄ and dq̄
    #q profile is set specifically for island coordinates.
    q, dq = q_prof(κ)

    #B = q̄ ∇κ×∇ᾱ + 2A ∇κ×∇τ
    
    B.B[1] = 0
    B.B[2] = 2*isl.A/ met.J[1]

    B.B[3] = q / met.J[1] 


    B.dB[2, 1] = - met.dJ[1] * 2 * isl.A / met.J[1]^2
    B.dB[2, 2] = - met.dJ[2] * 2 * isl.A / met.J[1]^2


    B.dB[3, 1] = dq / met.J[1] - q * met.dJ[1] / met.J[1]^2 
    B.dB[3, 2] = - met.dJ[2] * q / (met.J[1]^2)

    magnitude_B!(B, met)
    
    B.b[1] = B.B[1]/B.mag_B[1]
    B.b[2] = B.B[2]/B.mag_B[1]
    B.b[3] = B.B[3]/B.mag_B[1]

end

