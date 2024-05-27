
#making I better!

#start with the undampped case.
function new_compute_I!(I::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, met::MetT, B::BFieldT, n::Float64, δ::Float64)

    #without damping, I is a 3x3 matrix.
    #I = zeros(ComplexF64, 3, 3)

    #H is a length 9 array that is pre-allocated which is just a temp array.
    

    #need this for damping!
    H = zeros(ComplexF64, 9)

    for j in 1:3

        for i in 1:3
            H[j] += 1/met.J * (met.dJ[i] * (met.gu[i, j] - B.b[i]*B.b[j])
                + met.J * (met.dgu[i, j, i] - B.b[i]*B.db[j, i] - B.db[i, i]*B.b[j]))
        end
    end

    #Φ_ss
    H[4] = met.gu[1, 1]-B.b[1]*B.b[1]
    #Φ_sθ
    H[5] = met.gu[1, 2] - B.b[1]*B.b[2] + met.gu[2, 1] - B.b[2]*B.b[1] #symmetric could just double
    #Φ_sζ
    H[6] = met.gu[1, 3] - B.b[1]*B.b[3] + met.gu[3, 1] - B.b[3]*B.b[1]
    #Φ_θθ
    H[7] = met.gu[2, 2] - B.b[2]*B.b[2]
    #Φ_θζ
    H[8] = met.gu[2, 3] - B.b[2]*B.b[3] + met.gu[3, 2] - B.b[3]*B.b[2]
    #Φ_ζζ
    H[9] = met.gu[3, 3] - B.b[3]*B.b[3] 

    #display(met.J)
    #display(δ)

    I[:, :] = - (n * δ * met.J / B.mag_B^2 * 1.0im) .* H .* H' 

    #display(I)
    #I = -1.0 .* I
    #I .*= - n * δ * met.J / B.mag_B^2 * 1im
    """
    for j in 1:9, i in 1:9
        #for some reason the order of operation here seems to matter...
        #think this doesn't actually matter for us though.
        #hopefully we dont need to use float128...
        #this is especially susceptible to having met.J or H[i], H[j] first. which doesn't make sense but oh well
        #I guess we need to ignore differences of e-13 or more?? Seems wild.
        #with delta a tyoical value, these diffs are negligible, but speak of a wider issue...
        I[i, j] = - met.J * n / B.mag_B^2 * 1.0im * H[i] * H[j] * δ
        #I[i, j] = -n*δ*1im*H[i] * H[j]  * met.J / B.mag_B^2
    end
    """

    #I =  H .* H' #.*  (n * δ * met.J / B.mag_B^2)

    #I .*= n * δ * met.J / B.mag_B^2

    #this way around we are now getting ~e-22 errors in the real component. 
    #unclear if this is just going to be the problem with this code in general
    #Note, CKA, based on thesis, uses double precision.
    #this difference gives indistinguishable differences in computed frequency, at least for axel case.
    I[1:3, 1:3] += (met.J * n / B.mag_B^2) .* (met.gu[:, :] - B.b * B.b')  

    #do this second so that all values are changed first.
    #for i in 1:3, j in 1:3
    #    I[i, j] += met.J * n * (met.gu[i, j] - B.b[i] * B.b[j]) / B.mag_B^2
    #end

    return H
end