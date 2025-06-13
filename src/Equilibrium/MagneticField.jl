
"""
Struct for storing the magnetic field and related variables at a given coordinate.

### Fields
- B::Array{Float64} - Vector storing the magnetic field.
- b::Array{Float64} - Vector storing the normalised magnetic field.
- dB::Array{Float64, 2} V- ector storing the derivative of the magnetic field, second index refers to derivative coordinate.
- db::Array{Float64, 2} - Vector storing the derivative of the normalised magnetic field, second index refers to derivative coordinate.
- mag_B::Array{Float64} - Magnitude of the magnetic field. Stored as an array so struct is immutable.
- dmag_B::Array{Float64} - Derivative of the magnitude of B, index refers to derivative coordinate.
"""
struct BFieldT
    B :: Array{Float64, 1} 
    b :: Array{Float64, 1} 
    dB :: Array{Float64, 2} 
    db :: Array{Float64, 2} 
    mag_B :: Array{Float64, 1}
    dmag_B :: Array{Float64, 1} 
    function BFieldT()
        new(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), zeros(1), zeros(3))
    end
end



"""
    compute_B!(B::BFieldT, met::MetT, q_prof::Function, isl::IslandT, x1::Float64, x2::Float64, x3::Float64)

Function the fills the BFieldT struct based on an q-profile and and island for a given coordinate.
B is assumed to be in the form (B^r, B^x2, B^x3) with 
 - B^r = [A1(r) sin(m1_0 x2 + n1_0 x3) + A2(r) sin(m2_0 x2 + n2_0 x3)] / J
 - B^x2 = r/(J q) 
 - B^x3 = r / J

This version does not take a quadratic form, this means the behaviour near r=0 will be cooked, requires a restricted grid.
"""
function rad_compute_B!(B::BFieldT, met::MetT, q_prof::Function, isls::Array{IslandT}, x1::Float64, x2::Float64, x3::Float64)

    q, dq = q_prof(x1)


    #I think it is fine to just modify the B^r component, but we may want to confirm.
    #perhaps d amp is not a great name for the derivative.
    amp, damp = island_amplitude(x1)
    #amplitude seems to be a problem that will need to be fixed another time.
    #amp = 1.0
    #damp = 0.0 

    B.B[1] = 0.0
    B.dB[1, :] .= 0.0
    for isl in isls
        B.B[1] += isl.A * amp * isl.m0 * sin(isl.m0 * x2 + isl.n0 * x3)
        B.dB[1, 1] += isl.A * damp * isl.m0 * sin(isl.m0 * x2 + isl.n0 * x3)
        B.dB[1, 2] += isl.A * amp * isl.m0^2 * cos(isl.m0 * x2 + isl.n0 * x3)
        B.dB[1, 3] += isl.A * amp * isl.m0 * isl.n0 * cos(isl.m0 * x2 + isl.n0 * x3)
    end

    B.B[1] /= met.J[1]

    B.dB[1, 1] /= met.J[1]
    B.dB[1, 1] -= B.B[1] * met.dJ[1] / met.J[1]

    B.dB[1, 2] /= met.J[1]
    B.dB[1, 2] -= B.B[1] * met.dJ[2] / met.J[1]

    B.dB[1, 3] /= met.J[1]
    B.dB[1, 3] -= B.B[1] * met.dJ[3] / met.J[1] #typically dJ/dx3 = 0.

    #assumes B0=1
    #B.B[1] = 1 / (met.J[1]) * (isl.A * amp * isl.m0 * sin(arg) + isl2.A * amp * isl2.m0 * sin(arg2))
                
    B.B[2] = x1 / (met.J[1] * q) 
    B.B[3] = x1 / (met.J[1])

    #ignoring amp/damp for now
    #added it back in now!
    #B.dB[1, 1] = - met.dJ[1] / met.J[1]^2 * (isl.A * amp * isl.m0 * sin(arg) + isl2.A * amp * isl2.m0 * sin(arg2)) + 1 / (met.J[1]) * (isl.A * damp * isl.m0 * sin(arg) + isl2.A * damp  * isl2.m0 * sin(arg2))

    #B.dB[1, 2] = - met.dJ[2] / met.J[1]^2 * (isl.A * amp * isl.m0 * sin(arg) + isl2.A * amp * isl2.m0 * sin(arg2)) + 1 / met.J[1] * (isl.A * amp * isl.m0^2 * cos(arg) + isl2.A * amp * isl2.m0^2 * cos(arg2))

    #B.dB[1, 3] = + 1/ met.J[1] * (isl.A * amp * isl.m0 * isl.n0 * cos(arg) + isl2.A * amp * isl2.m0 * isl2.n0 * cos(arg2))

    #B.dB[1, 1] = ( - isl.A * amp * isl.m0 * sin(arg) * met.dJ[1] / met.J[1]^2
    #                + 1 / (met.J[1]) * isl.A * damp * isl.m0 * sin(arg))

    #B.dB[1, 2] = (1 / (met.J[1]) * isl.A * amp * isl.m0^2 * cos(arg)
    #                - isl.A * amp * isl.m0 * sin(arg) * met.dJ[2] / met.J[1]^2)

    #B.dB[1, 3] = isl.n0 / (met.J[1]) * isl.A * amp * isl.m0 * cos(arg)



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
    
    #unsure if there should be extra approximation here as Br is a pert and therefore small?
    B.b[1] = B.B[1]/B.mag_B[1]
    B.b[2] = B.B[2]/B.mag_B[1]
    B.b[3] = B.B[3]/B.mag_B[1]
end


"""
    magnitude_B!(B::BFieldT, met::MetT)

Function that computes the magnitude of B and its derivative once the other fields of BFieldT have been computed.
"""
function magnitude_B!(B::BFieldT, met::MetT)

    #B2 = 0.0
    #dB2 = zeros(Float64, 3)
    #changing this, hopefully not a problemo!
    B.mag_B .= 0.0
    B.dmag_B .= 0.0
    
    for i in 1:3, j in 1:3

        #note that this is actual |B|^2
        #B2 += B.B[i] * B.B[j] * met.gl[i, j]
        B.mag_B[1] += B.B[i] * B.B[j] * met.gl[i, j]

        for k in 1:3
            #should this be views? May need to look into that!
            #dB2[k] += (B.B[i] * B.B[j] * met.dgl[i, j, k] +
            #            B.B[i] * met.gl[i, j] * B.dB[j, k] +
            #            B.B[j] * met.gl[i, j] * B.dB[i, k])
            #needs to be modified to actually be the correct expression
            #We are just reusing memory
            B.dmag_B[k] += (B.B[i] * B.B[j] * met.dgl[i, j, k] +
                        B.B[i] * met.gl[i, j] * B.dB[j, k] +
                        B.B[j] * met.gl[i, j] * B.dB[i, k])
        end
    end

    B.mag_B[1] = sqrt(B.mag_B[1])
    B.dmag_B[:] = @. 1/(2*B.mag_B[1]) * B.dmag_B
    #B.mag_B = sqrt(B2)
    #B.dmag_B[:] = @. 1/(2*sqrt(B2)) * dB2

    for i in 1:3, j in 1:3
        B.db[j, i] = B.dB[j, i]/B.mag_B[1] - B.B[j]*B.dmag_B[i]/B.mag_B[1]^2
    end
end



"""

Function to determine the island amplitude. Needs work
"""
function island_amplitude(x1::Float64)

    #doesn't actually use the island at the moment, but probably will one day.
    #currently will only work for islands at r0=0.5.

    #ideally, this would be an input like q or density tbh

    #returns value and derivative.
    #return 4*x1*(1-x1), 4 - 8 * x1

    #case where this function is not included
    #useful to distinguish problemo's between GAM and between axis.
    #return 1, 0

    #we may want to try this for a flatter profile over the island
    return 1-16*(x1-1/2)^4, -64*(x1-1/2)^3

end


"""
    compute_B_isl!(B::BFieldT, met::MetT, isl::IslandT, κ::Float64, ᾱ::Float64, φ::Float64)

Computes the magnetic field for the island case. This requires a specific q-profile and a conversion between r and flux ψ.
"""
function compute_B_isl!(B::BFieldT, met::MetT, q_prof::Function, isl::IslandT, κ::Float64, ᾱ::Float64, τ::Float64)


    #I think to remove this function, we would have to swap from κ to either r̄ or χ
    #both are v annoying, these would also require a weird new q-profile.
    #and would make the metric extra cooked.
    #I think the way we do B is just a bit cooked, i.e. the normal case assumes Bθ has B_0 r or whatever.
    #think we have been assuming compute_B is a locked in function that can't be any different
    #but there is many a reason why this would want to change

    #should probably be q̄ and dq̄
    #note that this q profile is the same as toroidal case, but appears different in the magnetic field due to coordinates
    #this is the same as q appearing in other component if poloidal not toroidal flux is chosen as radial coordinate.
    q, dq = q_prof(κ)
    #K, E = Elliptic.ellipke(κ)

    #wot even is this function???
    #holy moly how did this ever work even a little bit.
    #this is from our case, this will not match Axel, but jopefully close enough...
    #A = 0.00015625000000000003
    #w = 0.05
    #A = 5.625e-5
    #w = 0.03
    #m0 = 2
    #n0 = -1
    #q = -w/(2*A*π*m0) * Elliptic.K(κ)
    #in built q-profile for island coordinates.
    #q = -isl.w/(2*isl.A*π*isl.m0) * K

    #this was a key peice!!!
    #this is just q for fuck sake, difference is the qprofile is at a different place!
    #maybe calling it the q-profile is perhaps a bit misleading then.
    #I guess this is the difference of using toroidal vs poloidal flux.
    #dψ̄dκ = isl.w * K / (isl.m0*π)
    #d2ψ̄dκ2 = isl.w / (isl.m0*π) * (E - (1-κ)*K) / (2*(1-κ)*κ)

    #dq = -w/(2*A*π*m0) * (Elliptic.E(κ) - (1-κ)*Elliptic.K(κ)) / (2*(1-κ)*κ)
    #dq = 0
    #dq = -isl.w/(2*isl.A*π*isl.m0) * (E - (1-κ)*K) / (2*(1-κ)*κ)

    #B = q̄ ∇κ×∇ᾱ + 2A ∇κ×∇τ
    #awkward units give this expression
    #but this is mainly so we can define our metric in terms of κ, instead of r̄ or χ, both of which would have nicer B expressions.
    #much simpler expressions now we understand our damn coordinates
    #it may be worth considering what the metric would look like in χ or r̄
    #would be a disaster nevermind.
    
    B.B[1] = 0
    B.B[2] = 2*isl.A/ met.J[1]

    B.B[3] = q / met.J[1] 


    B.dB[2, 1] = - met.dJ[1] * 2 * isl.A / met.J[1]^2
    B.dB[2, 2] = - met.dJ[2] * 2 * isl.A / met.J[1]^2


    B.db[3, 1] = dq / met.J[1] - q * met.dJ[1] / met.J[1]^2
    B.dB[3, 2] = - met.dJ[2] * q / (met.J[1]^2)


    magnitude_B!(B, met)
    
    B.b[1] = B.B[1]/B.mag_B[1]
    B.b[2] = B.B[2]/B.mag_B[1]
    B.b[3] = B.B[3]/B.mag_B[1]




end




"""
    compute_B!(B::BFieldT, met::MetT, q_prof::Function, isl::IslandT, x1::Float64, x2::Float64, x3::Float64)

Function the fills the BFieldT struct based on an q-profile and and island for a given coordinate.
B is assumed to be in the form (B^r, B^x2, B^x3) with 
 - B^r = A m_0 r^2 (1-r) sin(m_0 x2 - n_0 x3) / J
 - B^x2 = r/(J q) + A (1-2r) cos(m_0 x2 - n_0 x3) / J
 - B^x3 = r / J

This version does not take a quadratic form, this means the behaviour near r=0 will be cooked, requires a restricted grid.
This has been changed to ψ, needs verification.
"""
function flux_compute_B!(B::BFieldT, met::MetT, q_prof::Function, isl::IslandT, ψ::Float64, x2::Float64, x3::Float64)

    q, dq = q_prof(ψ)

    arg = isl.m0 * x2 + isl.n0 * x3

    #assumes B0=1
    B.B[1] = 1 / (met.J[1]) * isl.A * isl.m0 * sin(arg)
                
    B.B[2] = 1 / (met.J[1] * q) 
    B.B[3] = 1 / (met.J[1])

    B.dB[1, 1] = ( - isl.A * isl.m0 * sin(arg) * met.dJ[1] / met.J[1]^2)

    B.dB[1, 2] = (1 / (met.J[1]) * isl.A * isl.m0^2 * cos(arg)
                    - isl.A * isl.m0 * sin(arg) * met.dJ[2] / met.J[1]^2)

    B.dB[1, 3] = isl.n0 / (met.J[1]) * isl.A * isl.m0 * cos(arg)



    B.dB[2, 1] = ( - (1 / q) * met.dJ[1] / met.J[1]^2
                    - 1 * dq /(met.J[1] * q^2) )
                    
    B.dB[2, 2] = - (1 / q ) * met.dJ[2] / met.J[1]^2
    


    B.dB[3, 1] = ( - met.dJ[1]) / met.J[1]^2
    B.dB[3, 2] = - met.dJ[2] / met.J[1]^2
    

    #note this also does B.db
    magnitude_B!(B, met)
    
    #unsure if there should be extra approximation here as Br is a pert and therefore small?
    B.b[1] = B.B[1]/B.mag_B
    B.b[2] = B.B[2]/B.mag_B
    B.b[3] = B.B[3]/B.mag_B
end


