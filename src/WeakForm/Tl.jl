
"""
    function Tl!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, met::MetT, B::BFieldT, C::Array{Float64, 2}, D::Array{Float64, 2}, T::Array{Float64, 2})

Computes the laplacian like term of W.
"""
function Tl!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, met::MetT, B::BFieldT, C::Array{Float64, 2}, D::Array{Float64, 2}, T::Array{Float64, 2})

    compute_C!(B, C)
    #compute_C_old!(B, C)

    mul!(T, C', D)

    mul!(W, T, C)

    W .*= met.J * B.mag_B^2
    #option for old C
    #W .*= met.J #* B.mag_B^2
end


"""
    function compute_C!(B::BFieldT, C::Array{Float64, 2})

Computes the C matrix used to compute Tl.
"""
function compute_C!(B::BFieldT, C::Array{Float64, 2})

    
    #C_i^j = ∂_i (B^j/B^2)
    #for the first 3x3 block.
    for i in 1:3
        for j in 1:3

            #could be very wrong
            C[j, i] = B.dB[i, j] / B.mag_B^2 - 2 * B.B[i] * B.dmag_B[j] / B.mag_B^3 

        end
    end

    #Rest of C is given by B^i/B^2
    #weird pattern comes from only represneting double derivatives 
    #once.

    C[1, 4] = B.B[1] / B.mag_B^2

    C[1, 5] = B.B[2] / B.mag_B^2
    C[2, 5] = B.B[1] / B.mag_B^2

    C[1, 6] = B.B[3] / B.mag_B^2
    C[3, 6] = B.B[1] / B.mag_B^2

    C[2, 7] = B.B[2] / B.mag_B^2

    C[2, 8] = B.B[3] / B.mag_B^2
    C[3, 8] = B.B[2] / B.mag_B^2

    C[3, 9] = B.B[3] / B.mag_B^2


end

