
"""
    function compute_W!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, met::MetT, B::BFieldT, n::Float64, ωcap2::Float64, tm::TM)

Computes the W matrix for the weak form at a single coordinate.
"""
function compute_W!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, met::MetT, B::BFieldT, n::Float64, ωcap2::Float64, tm::TM)


   
    #compute the laplacian like term
    Tl!(W, met, B, tm.C, tm.D, tm.T)
    
    #display("W is pre Tj")
    #display(W[1:3, 1:3])
    #display(W[1:3, 4:6])
    #display(W[1:3, 7:9])

    #compute the current term.
    Tj!(W, met, B, tm.Γ, tm.dΓ, tm.K)
    


    #display("W is")
    #display(W[1:3, 1:3])
    #display(W[1:3, 4:6])
    #display(W[1:3, 7:9])
    #TODO
    #work in progress.
    #W[1:3, 1:3] += tm.D .* (ωcap2 * n *  met.J / B.mag_B^2)

end


"""
    function compute_isl_W!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, met::MetT, B::BFieldT, tm::TM)

Computes the W matrix in the case of island coordinates. In this case, the current term is excluded.
"""
function compute_isl_W!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, met::MetT, B::BFieldT, tm::TM)


    Tl!(W, met, B, tm.C, tm.D, tm.T)
    
end




#######################################
#old stuff.

#for comparison...
function compute_C_old!(B::BFieldT, C::Array{Float64, 2})

    #if we change this, perhaps we won't need to store \bb anymore?
    #although I think it is used in current term.

    for i in 1:3
        for j in 1:3

            #could be very wrong
            C[j, i] = B.db[i, j]

        end
    end

    C[1, 4] = B.b[1]

    C[1, 5] = B.b[2]
    C[2, 5] = B.b[1]

    C[1, 6] = B.b[3]
    C[3, 6] = B.b[1]

    C[2, 7] = B.b[2]

    C[2, 8] = B.b[3]
    C[3, 8] = B.b[2]

    C[3, 9] = B.b[3]


end

#think this is actually not needed
#as island case does not compute dB or dg, so current term will always be zero regardless.
#also, we may eventually want those additions,
#however, this should be a bit faster.



#appears to be giving the same results as other cases.
#identical results to earlier
#but we have changed compute B and are now getting different results
"""
    compute_Tj(met::MetT, B::BFieldT)

Computes the current term contribution to W.
"""
function compute_Tj(met::MetT, B::BFieldT)
    Tj = zeros(9, 9) #9x9 is overkill, as Tj[>3, >3] will always be zero.

    Γ = zeros(Float64, 3, 3)
    dΓ = zeros(Float64, 3, 3, 3) #last ind is deriv.

    #Γ_i^j = δ_i^j - g_{im}b^m b^j
    for i in 1:3, j in 1:3
        Γ[i, j] = dδ(i, j) - dot(met.gl[i, :], B.b[:])*B.b[j]
    end

    for i in 1:3, j in 1:3, k in 1:3

        dΓ[i, j, k] = (- B.b[j] * dot(B.b[:], met.dgl[i, :, k]) 
                        - dot(met.gl[i, :], B.b[:]) * B.db[j, k]
                        - B.b[j] * dot(met.gl[i, :], B.db[:, k]))
    end
    #Tj  is given by the expression
    #Γ_i^n ∂_nΨ 1/J ϵ^{ijk}(Γ_k^q∂_j∂_qΦ + ∂_j(Γ_k^q)∂_qΦ) + Γ_i^n ∂_nΦ 1/J ϵ^{ijk}(Γ_k^q∂_j∂_qΨ + ∂_j(Γ_k^q)∂_qΨ)

    #first we consdier the 3x3 of Tj, which only has first derivative terms.
    #Γ_i^n ∂_nΨ 1/J ϵ^{ijk}∂_j(Γ_k^q)∂_qΦ + Γ_i^n ∂_nΦ 1/J ϵ^{ijk}∂_j(Γ_k^q)∂_qΨ
    for n in 1:3, q in 1:3
        val = 0
        for i in 1:3, j in 1:3, k in 1:3
            val += Γ[i, n] / met.J * lct[i, j, k] * dΓ[k, q, j]
        end
        Tj[n, q] += val
        Tj[q, n] += val
    end

    #now we consider the second derivates.

    for n in 1:3

        #this still might be improvable, but is much more concise.
        #why the fk is this called dds.
        dds = zeros(6)

        for i in 1:3, k in 1:3
            dds[1] += Γ[i, n] * lct[i, 1, k] * Γ[k, 1] / met.J 

            dds[2] += Γ[i, n] * (lct[i, 1, k] * Γ[k, 2] + lct[i, 2, k] * Γ[k, 1]) / met.J

            dds[3] += Γ[i, n] * (lct[i, 1, k] * Γ[k, 3] + lct[i, 3, k] * Γ[k, 1]) / met.J

            dds[4] += Γ[i, n] * lct[i, 2, k] * Γ[k, 2] / met.J

            dds[5] += Γ[i, n] * (lct[i, 2, k] * Γ[k, 3] + lct[i, 3, k] * Γ[k, 2]) / met.J

            dds[6] += Γ[i, n] * lct[i, 3, k] * Γ[k, 3] / met.J
        end

        Tj[n, 4:9] .+= dds

        Tj[4:9, n] .+= dds
    end
         

    return Tj

end


#no longer used.
function compute_Tl(met::MetT, B::BFieldT)

    #love that we are still doing this...
    C = zeros(3, 9)
    D = zeros(3, 3) #maybe shouldn't recompute this.
    Tl = zeros(9, 9)

    #if we change this, perhaps we won't need to store \bb anymore?
    #although I think it is used in current term.

    for i in 1:3
        for j in 1:3

            #could be very wrong
            C[j, i] = B.dB[i, j] / B.mag_B^2 - 2 * B.B[i] * B.dmag_B[j] / B.mag_B^3 

            D[i, j] = met.gu[i, j] - B.b[i]*B.b[j]


        end
    end

    C[1, 4] = B.B[1] / B.mag_B^2

    C[1, 5] = B.B[2] / B.mag_B^2
    C[2, 5] = B.B[1] / B.mag_B^2

    C[1, 6] = B.B[3] / B.mag_B^2
    C[3, 6] = B.B[1] / B.mag_B^2

    C[2, 7] = B.B[2] / B.mag_B^2

    C[2, 8] = B.B[3] / B.mag_B^2
    C[3, 8] = B.B[2] / B.mag_B^2

    C[3, 9] = B.B[3] / B.mag_B^2



    for μ in 1:9

        for ν in 1:9

            for i in 1:3
                for j in 1:3

                    #would this be faster as a matrix multiplication?
                    #maybe more memory efficient?
                    #then we don't need to set Tl to zeros??
                    Tl[μ, ν] += C[i, μ] * D[i, j] * C[j, ν]
                end
            end
        end
    end

    return Tl


end



#this uses the weak form version in CKA thesis, this assumes the derivative of |B| is small compared to other things, which is not really a valid assumption for our case.
function compute_Tl_old(met::MetT, B::BFieldT)

    C = zeros(3, 9)
    D = zeros(3, 3) #maybe shouldn't recompute this.
    Tl = zeros(9, 9)

    for i in 1:3
        for j in 1:3
            C[j, i] = B.db[i, j] 

            D[i, j] = met.gu[i, j] - B.b[i]*B.b[j]
            
        end
    end
    """
    #not sure this will work as we want, maybe this is the problem with other term?
    for i in 1:3
        for j in 1:3
            C[i, d_to_ind(i, j)] += B.b[j]
        end
    end

    #can probably do this automatically with our d_to_ind func
    """
    C[1, 4] = B.b[1]

    C[1, 5] = B.b[2]
    C[2, 5] = B.b[1]

    C[1, 6] = B.b[3]
    C[3, 6] = B.b[1]

    C[2, 7] = B.b[2]

    C[2, 8] = B.b[3]
    C[3, 8] = B.b[2]

    C[3, 9] = B.b[3]


    for μ in 1:9

        for ν in 1:9

            for i in 1:3
                for j in 1:3

                    Tl[μ, ν] += C[i, μ] * D[i, j] * C[j, ν]
                end
            end
        end
    end

    return Tl

end
