
"""
Struct for storing the magnetic field and related variables at a given coordinate.

### Fields
- B::Array{Float64} - Vector storing the magnetic field.
- b::Array{Float64} - Vector storing the normalised magnetic field.
- dB::Array{Float64, 2} V- ector storing the derivative of the magnetic field, second index refers to derivative coordinate.
- db::Array{Float64, 2} - Vector storing the derivative of the normalised magnetic field, second index refers to derivative coordinate.
- mag_B::Float64 - Magnitude of the magnetic field.
- dmag_B::Array{Float64} - Derivative of the magnitude of B, index refers to derivative coordinate.
"""
mutable struct BFieldT
    B :: Array{Float64, 1} 
    b :: Array{Float64, 1} 
    dB :: Array{Float64, 2} 
    db :: Array{Float64, 2} 
    mag_B :: Float64 
    dmag_B :: Array{Float64, 1} 
    function BFieldT()
        new(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), 0.0, zeros(3))
    end
end


"""
    init_empty_B()

Initialises empty BFieldT struct.
"""
function init_empty_B()

    return BFieldT(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), 0.0, zeros(3))
end


"""
    compute_B!(B::BFieldT, met::MetT, q_prof::Function, isl::IslandT, r::Float64, θ::Float64, ζ::Float64)

Function the fills the BFieldT struct based on an q-profile and and island for a given coordinate.
B is assumed to be in the form (B^r, B^θ, B^ζ) with 
 - B^r = A m_0 r^2 (1-r) sin(m_0 θ - n_0 ζ) / J
 - B^θ = r/(J q) + A (1-2r) cos(m_0 θ - n_0 ζ) / J
 - B^ζ = r / J

This version does not take a quadratic form, this means the behaviour near r=0 will be cooked, requires a restricted grid.
"""
function compute_B!(B::BFieldT, met::MetT, q_prof::Function, isl::IslandT, isl2::IslandT, r::Float64, θ::Float64, ζ::Float64)

    q, dq = q_prof(r)

    arg = isl.m0 * θ + isl.n0 * ζ

    arg2 = isl2.m0 * θ + isl2.n0 * ζ


    #I think it is fine to just modify the B^r component, but we may want to confirm.
    #perhaps d amp is not a great name for the derivative.
    amp, damp = island_amplitude(r, isl)
    #amplitude seems to be a problem that will need to be fixed another time.
    #amp = 1.0
    #damp = 0.0 

    #assumes B0=1
    B.B[1] = 1 / (met.J) * (isl.A * amp * isl.m0 * sin(arg) + isl2.A * amp * isl2.m0 * sin(arg2))
                
    B.B[2] = r / (met.J * q) 
    B.B[3] = r / (met.J)

    #ignoring amp/damp for now
    #added it back in now!
    B.dB[1, 1] = - met.dJ[1] / met.J^2 * (isl.A * amp * isl.m0 * sin(arg) + isl2.A * amp * isl2.m0 * sin(arg2)) + 1 / (met.J) * (isl.A * damp * isl.m0 * sin(arg) + isl2.A * damp  * isl2.m0 * sin(arg2))

    B.dB[1, 2] = - met.dJ[2] / met.J^2 * (isl.A * amp * isl.m0 * sin(arg) + isl2.A * amp * isl2.m0 * sin(arg2)) + 1 / met.J * (isl.A * amp * isl.m0^2 * cos(arg) + isl2.A * amp * isl2.m0^2 * cos(arg2))

    B.dB[1, 3] = + 1/ met.J * (isl.A * amp * isl.m0 * isl.n0 * cos(arg) + isl2.A * amp * isl2.m0 * isl2.n0 * cos(arg2))

    #B.dB[1, 1] = ( - isl.A * amp * isl.m0 * sin(arg) * met.dJ[1] / met.J^2
    #                + 1 / (met.J) * isl.A * damp * isl.m0 * sin(arg))

    #B.dB[1, 2] = (1 / (met.J) * isl.A * amp * isl.m0^2 * cos(arg)
    #                - isl.A * amp * isl.m0 * sin(arg) * met.dJ[2] / met.J^2)

    #B.dB[1, 3] = isl.n0 / (met.J) * isl.A * amp * isl.m0 * cos(arg)



    B.dB[2, 1] = (1 / (met.J * q) 
                    - (r / q) * met.dJ[1] / met.J^2
                    - r * dq /(met.J * q^2) )
                    
    B.dB[2, 2] = - (r / q ) * met.dJ[2] / met.J^2
    


    B.dB[3, 1] = (met.J - r * met.dJ[1]) / met.J^2
    B.dB[3, 2] = -r*met.dJ[2]/met.J^2
    

    #note this also does B.db
    magnitude_B!(B, met)
    
    #unsure if there should be extra approximation here as Br is a pert and therefore small?
    B.b[1] = B.B[1]/B.mag_B
    B.b[2] = B.B[2]/B.mag_B
    B.b[3] = B.B[3]/B.mag_B
end


function island_amplitude(r::Float64, isl::IslandT)

    #doesn't actually use the island at the moment, but probably will one day.
    #currently will only work for islands at r0=0.5.

    #ideally, this would be an input like q or density tbh

    #returns value and derivative.
    #return 4*r*(1-r), 4 - 8 * r

    #case where this function is not included
    #useful to distinguish problemo's between GAM and between axis.
    return 1, 0

    if r > 0.3 #this will assume the island is at r=0.5
        return 1, 0
    elseif r > 0.2
        #function chosen to be zero at r=0.2, 1 at r=0.3, and gradient of 0 at r=0.3.
        #subject to change obvs.
        return -100*r^2 + 60*r - 8, -200*r + 60
    else

        return 0, 0
    end

end

#simplified function that only computes the contravarient components
#useful for QFM surfaces.
function compute_B!(B::Array{Float64}, J::Float64, q_prof::Function, isl::IslandT, isl2::IslandT, r::Float64, θ::Float64, ζ::Float64)
    
    q, _ = q_prof(r)

    arg = isl.m0 * θ + isl.n0 * ζ

    arg2 = isl2.m0 * θ + isl2.n0 * ζ

    #assumes B0=1
    B[1] = 1 / J * (isl.A * isl.m0 * sin(arg) + isl2.A  * isl2.m0 * sin(arg2))
                
    B[2] = r / (J * q) 
    B[3] = r / (J)

end

using Elliptic

#computes B when in island coordinates. very different who knew.
function compute_B_isl!(B::BFieldT, met::MetT, isl::IslandT, κ::Float64, ᾱ::Float64, φ::Float64)

    #q, dq = q_prof(r)


    #this is from our case, this will not match Axel, but jopefully close enough...
    #A = 0.00015625000000000003
    #w = 0.05
    #A = 5.625e-5
    #w = 0.03
    #m0 = 2
    #n0 = -1
    #q = -w/(2*A*π*m0) * Elliptic.K(κ)
    #in built q-profile for island coordinates.
    q = -isl.w/(2*isl.A*π*isl.m0) * Elliptic.K(κ)

    #this was a key peice!!!
    dψ̄dκ = isl.w * Elliptic.K(κ) / (isl.m0*π)

    #dq = -w/(2*A*π*m0) * (Elliptic.E(κ) - (1-κ)*Elliptic.K(κ)) / (2*(1-κ)*κ)
    dq = 0
    
    B.B[1] = 0
    #Axel just has R0 here, which doesn't make anysense!
    B.B[2] = 1/(q * met.J) * dψ̄dκ

    B.B[3] = 1 / met.J * dψ̄dκ


    B.dB[2, 1] = -dq / (q^2*met.J) - met.dJ[1] / (q*met.J^2)
    B.db[2, 2] = - met.dJ[2] / (q*met.J^2)


    B.dB[3, 1] = - met.dJ[1] / (met.J^2)
    B.dB[3, 2] = - met.dJ[2] / (met.J^2)


    magnitude_B!(B, met)
    
    B.b[1] = B.B[1]/B.mag_B
    B.b[2] = B.B[2]/B.mag_B
    B.b[3] = B.B[3]/B.mag_B




end




"""
    compute_B!(B::BFieldT, met::MetT, q_prof::Function, isl::IslandT, r::Float64, θ::Float64, ζ::Float64)

Function the fills the BFieldT struct based on an q-profile and and island for a given coordinate.
B is assumed to be in the form (B^r, B^θ, B^ζ) with 
 - B^r = A m_0 r^2 (1-r) sin(m_0 θ - n_0 ζ) / J
 - B^θ = r/(J q) + A (1-2r) cos(m_0 θ - n_0 ζ) / J
 - B^ζ = r / J

This version does not take a quadratic form, this means the behaviour near r=0 will be cooked, requires a restricted grid.
This has been changed to ψ, needs verification.
"""
function flux_compute_B!(B::BFieldT, met::MetT, q_prof::Function, isl::IslandT, ψ::Float64, θ::Float64, ζ::Float64)

    q, dq = q_prof(ψ)

    arg = isl.m0 * θ + isl.n0 * ζ

    #assumes B0=1
    B.B[1] = 1 / (met.J) * isl.A * isl.m0 * sin(arg)
                
    B.B[2] = 1 / (met.J * q) 
    B.B[3] = 1 / (met.J)

    B.dB[1, 1] = ( - isl.A * isl.m0 * sin(arg) * met.dJ[1] / met.J^2)

    B.dB[1, 2] = (1 / (met.J) * isl.A * isl.m0^2 * cos(arg)
                    - isl.A * isl.m0 * sin(arg) * met.dJ[2] / met.J^2)

    B.dB[1, 3] = isl.n0 / (met.J) * isl.A * isl.m0 * cos(arg)



    B.dB[2, 1] = ( - (1 / q) * met.dJ[1] / met.J^2
                    - 1 * dq /(met.J * q^2) )
                    
    B.dB[2, 2] = - (1 / q ) * met.dJ[2] / met.J^2
    


    B.dB[3, 1] = ( - met.dJ[1]) / met.J^2
    B.dB[3, 2] = - met.dJ[2] / met.J^2
    

    #note this also does B.db
    magnitude_B!(B, met)
    
    #unsure if there should be extra approximation here as Br is a pert and therefore small?
    B.b[1] = B.B[1]/B.mag_B
    B.b[2] = B.B[2]/B.mag_B
    B.b[3] = B.B[3]/B.mag_B
end



"""
    compute_B!(B::BFieldT, met::MetT, q_prof::Function, isl::IslandT, r::Float64, θ::Float64, ζ::Float64)

Function the fills the BFieldT struct based on an q-profile and and island for a given coordinate.
B is assumed to be in the form (B^r, B^θ, B^ζ) with 
 - B^r = A m_0 r^2 (1-r) sin(m_0 θ - n_0 ζ) / J
 - B^θ = r/(J q) + A (1-2r) cos(m_0 θ - n_0 ζ) / J
 - B^ζ = r / J

This version does not take a quadratic form, this means the behaviour near r=0 will be cooked, requires a restricted grid.
"""
function constant_amp_compute_B!(B::BFieldT, met::MetT, q_prof::Function, isl::IslandT, r::Float64, θ::Float64, ζ::Float64)

    q, dq = q_prof(r)

    arg = isl.m0 * θ + isl.n0 * ζ

    #assumes B0=1
    B.B[1] = 1 / (met.J) * isl.A * isl.m0 * sin(arg)
                
    B.B[2] = r / (met.J * q) 
    B.B[3] = r / (met.J)

    B.dB[1, 1] = ( - isl.A * isl.m0 * sin(arg) * met.dJ[1] / met.J^2)

    B.dB[1, 2] = (1 / (met.J) * isl.A * isl.m0^2 * cos(arg)
                    - isl.A * isl.m0 * sin(arg) * met.dJ[2] / met.J^2)

    B.dB[1, 3] = isl.n0 / (met.J) * isl.A * isl.m0 * cos(arg)



    B.dB[2, 1] = (1 / (met.J * q) 
                    - (r / q) * met.dJ[1] / met.J^2
                    - r * dq /(met.J * q^2) )
                    
    B.dB[2, 2] = - (r / q ) * met.dJ[2] / met.J^2
    


    B.dB[3, 1] = (met.J - r * met.dJ[1]) / met.J^2
    B.dB[3, 2] = -r*met.dJ[2]/met.J^2
    

    #note this also does B.db
    magnitude_B!(B, met)
    
    #unsure if there should be extra approximation here as Br is a pert and therefore small?
    B.b[1] = B.B[1]/B.mag_B
    B.b[2] = B.B[2]/B.mag_B
    B.b[3] = B.B[3]/B.mag_B
end




"""
    compute_B!(B::BFieldT, met::MetT, q_prof::Function, isl::IslandT, r::Float64, θ::Float64, ζ::Float64)

Function the fills the BFieldT struct based on an q-profile and and island for a given coordinate.
B is assumed to be in the form (B^r, B^θ, B^ζ) with 
 - B^r = A m_0 r^2 (1-r) sin(m_0 θ - n_0 ζ) / J
 - B^θ = r/(J q) + A (1-2r) cos(m_0 θ - n_0 ζ) / J
 - B^ζ = r / J
"""
function compute_B_quadratic!(B::BFieldT, met::MetT, q_prof::Function, isl::IslandT, r::Float64, θ::Float64, ζ::Float64)

    q, dq = q_prof(r)

    arg = isl.m0 * θ + isl.n0 * ζ

    #assumes B0=1
    B.B[1] = 1 / (met.J) * isl.A * isl.m0 * r^2 * (1-r) * sin(arg)
                
    B.B[2] = r / (met.J * q) + isl.A / met.J * (1-2*r) * cos(arg)
    B.B[3] = r / (met.J)

    B.dB[1, 1] = (1 / (met.J) * isl.A * isl.m0 * (2*r-3*r^2) * sin(arg) 
                    - isl.A * isl.m0 * r^2 * (1-r) * sin(arg) * met.dJ[1] / met.J^2)

    B.dB[1, 2] = (1 / (met.J) * isl.A * isl.m0^2 * r^2 * (1-r) * cos(arg)
                    - isl.A * isl.m0 * r^2 * (1-r) * sin(arg) * met.dJ[2] / met.J^2)

    B.dB[1, 3] = isl.n0 / (met.J) * isl.A * isl.m0 * r^2 * (1-r) * cos(arg)



    B.dB[2, 1] = (1 / (met.J * q) - 2 * isl.A / met.J * cos(arg)
                    - (r / q - isl.A * (1-2*r) * cos(arg)) * met.dJ[1] / met.J^2
                    - r * dq /(met.J * q^2) )
                    
    B.dB[2, 2] = (- isl.m0 * isl.A / met.J * (1-2*r) * sin(arg)
                    - (r / q + isl.A * (1-2*r) * cos(arg)) * met.dJ[2] / met.J^2)
    
    #doing this correctly has changed a few things...
    B.dB[2, 3] = isl.n0 * isl.A / met.J * (1-2*r) * sin(arg)



    B.dB[3, 1] = (met.J - r * met.dJ[1]) / met.J^2
    B.dB[3, 2] = -r*met.dJ[2]/met.J^2
    

    #note this also does B.db
    magnitude_B!(B, met)
    
    #unsure if there should be extra approximation here as Br is a pert and therefore small?
    B.b[1] = B.B[1]/B.mag_B
    B.b[2] = B.B[2]/B.mag_B
    B.b[3] = B.B[3]/B.mag_B
end



"""
    magnitude_B!(B::BFieldT, met::MetT)

Function that computes the magnitude of B and its derivative once the other fields of BFieldT have been computed.
"""
function magnitude_B!(B::BFieldT, met::MetT)

    B2 = 0.0
    dB2 = zeros(Float64, 3)
    #
    for i in 1:3, j in 1:3

        B2 += B.B[i] * B.B[j] * met.gl[i, j]

        for k in 1:3
            #should this be views? May need to look into that!
            dB2[k] += (B.B[i] * B.B[j] * met.dgl[i, j, k] +
                        B.B[i] * met.gl[i, j] * B.dB[j, k] +
                        B.B[j] * met.gl[i, j] * B.dB[i, k])
        end
    end

    B.mag_B = sqrt(B2)
    B.dmag_B[:] = @. 1/(2*sqrt(B2)) * dB2

    for i in 1:3, j in 1:3
        B.db[j, i] = B.dB[j, i]/B.mag_B - B.B[j]*B.dmag_B[i]/B.mag_B^2
    end
end


"""
    compute_island_B!(B::BFieldT, met::MetT, isl::ContIslandT, ψ::Float64, α::Float64)

Function the fills the BFieldT struct specifically for the case of the island continuum.
This employs special island flux coordinates, uses a specific q-profile
and only requires the magnitude of B.
"""
function compute_island_B!(B::BFieldT, met::MetT, isl, ψ::Float64, α::Float64)

    #specifc q-profile required for the coordinate transformation.
    q = 1 / (1 / isl.q0 - isl.qp /isl.q0^2 * (ψ - isl.ψ0))

    #may want a factor or 1/4 in here to make this match our other cases.
    B.B[1] = isl.A * isl.m0 *sin(isl.m0 * α)
    B.B[2] = 1 / (met.J * q)
    B.B[3] = 1 / met.J 


    B2 = 0

    for i in 1:3, j in 1:3

        B2 += B.B[i] * met.gl[i, j] * B.B[j]
    end
    B.mag_B = sqrt(B2)

end
