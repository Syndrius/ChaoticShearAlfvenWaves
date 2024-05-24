
"""
Computes the I matrix for the weak form at a single coordinate. Still a bit unoptimised. See thesis for details on what is being computed.

# Args 
- I Stupid type due to using views.
- met::MetT Struct containing the metric information.
- B::BfieldT Struct containing the magnetic field information.
- n::Float64 Density.
- δ::Float64 Artifical resistivity, for damping calculations.
"""
function compute_I!(I::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, met::MetT, B::BFieldT, n::Float64, δ::Float64)

    #fill!(I, 0.0 + 0.0im)
    #display(I)

    #not sure how to deal with this yet!
    H = zeros(Float64, 9)
    for j=1:3, i=1:3
        #maybe this would be faster with built in matric operations??
        #I[i, j] = n*(met.gu[i, j] - B.b[i]*B.b[j])/B.mag_B^2 * met.J
        H[j] += 1/met.J * ( met.dJ[i]*(met.gu[i, j] - B.b[i]*B.b[j])
            + met.J*(met.dgu[i, j, i] - B.b[i]*B.db[j, i] - B.db[i, i] * B.b[j]))
    end
    #display(I)

    #something is wrong with this!!!! damn daniel
    H[4] = met.gu[1, 1]-B.b[1]*B.b[1]
    #ϕ_rθ, i=r, j=θ + i=θ, j=r
    H[5] = met.gu[1, 2]-B.b[1]*B.b[2] + met.gu[2, 1]-B.b[2]*B.b[1]
    #ϕ_rζ, i=r, j=ζ + i=ζ, j=r
    H[6] = met.gu[1, 3]-B.b[1]*B.b[3] + met.gu[3, 1]-B.b[3]*B.b[1]
    #ϕ_θθ, i=θ, j=θ
    H[7] = met.gu[2, 2]-B.b[2]*B.b[2] 
    #ϕ_θζ, i=θ, j=ζ + i=ζ, j=θ
    H[8] = met.gu[2, 3]-B.b[2]*B.b[3] + met.gu[3, 2]-B.b[3]*B.b[2]
    #ϕ_ζζ, i=ζ, j=ζ
    H[9] = met.gu[3, 3]-B.b[3]*B.b[3] 


    #probably want to do this as an outer product perhaps.
    for j=1:9, i=1:9
        I[i, j] = -n*δ*H[i] * H[j]  * met.J / B.mag_B^2*1.0im
    end
    #ugly fix!
    #so I is not being reset to zero for some reason
    #this forces all 9x9 entries to be filled by the damping term first, then we modify the un-damped case.
    for j=1:3, i=1:3
        I[i, j] += n*(met.gu[i, j] - B.b[i]*B.b[j]) * met.J /B.mag_B^2
    end

    return H
end