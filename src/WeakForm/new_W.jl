
function new_compute_W!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, met::MetT, B::BFieldT)

    #display("New guy")
    #start with the laplace like term
    W[:, :] = new_compute_Tl(met, B) .* met.J
    
    W[:, :] -= new_compute_Tj(met, B) .* met.J .* new_jparonB(met, B) ./ 2
    #Now getting identical results, at least in Axel case.
    #this implies it is not our implementation of this that is wrong
    #rather it is what we are implementing??
    #println(W)

end


#this looks to be identical.
function new_jparonB(met, B)
    #now you know we are in trouble!
    lc = zeros(3, 3, 3)
    lc[1, 2, 3] = 1
    lc[2, 3, 1] = 1
    lc[3, 1, 2] = 1
    lc[3, 2, 1] = -1
    lc[1, 3, 2] = -1
    lc[2, 1, 3] = -1

    J = zeros(3)
    for i in 1:3, j in 1:3, k in 1:3

        J[i] += 1 / met.J * lc[i, j, k] * (dot(met.gl[k, :], B.dB[:, j]) + dot(met.dgl[k, :, j], B.B[:]))
    end


    jpar = 0
    for i in 1:3, j in 1:3
        jpar += met.gl[i, j] * B.b[i] * J[j]
    end
    #display(jpar/B.mag_B)
    return jpar/B.mag_B

end

function dδ(i, j)
    if i==j
        return 1.0
    else
        return 0.0
    end
end


function new_compute_Tj(met, B)
    Tj = zeros(9, 9)

    #now you know we are in trouble!
    lc = zeros(3, 3, 3)
    lc[1, 2, 3] = 1.0
    lc[2, 3, 1] = 1.0
    lc[3, 1, 2] = 1.0
    lc[3, 2, 1] = -1.0
    lc[1, 3, 2] = -1.0
    lc[2, 1, 3] = -1.0

    #TODO!

    #lets just condier first case to see how cooked this will be...

    #maybe we consider just the 3x3 for now.

    Γ = zeros(3, 3)
    dΓ = zeros(3, 3, 3) #last ind is deriv.

    for i in 1:3, j in 1:3
        Γ[i, j] = dδ(i, j) - dot(met.gl[i, :], B.b[:])*B.b[j]
    end

    for i in 1:3, j in 1:3, k in 1:3

        dΓ[i, j, k] = (- B.b[j] * dot(B.b[:], met.dgl[i, :, k]) 
                        - dot(met.gl[i, :], B.b[:]) * B.db[j, k]
                        - B.b[j] * dot(met.gl[i, :], B.db[:, k]))
    end

    #Γ seems to be hte same in both methods. may or may not be good.
    #display(Γ)
    #display(dΓ)

    #defo some large amount of symmetry we ccould take advantage of here, start simple though!
    for i in 1:3, j in 1:3, k in 1:3
        #need to recheck these bad bois with both terms considered at once.
        Tj[1, 1] += Γ[i, 1] * lc[i, j, k] * dΓ[k, 1, j] / met.J * 2 #+ Γ[k, 1] * lc[i, j, k] * dΓ[i, 1, j] / met.J #from each side??
        Tj[1, 2] += Γ[i, 1] * lc[i, j, k] * dΓ[k, 2, j] / met.J + Γ[i, 2] * lc[i, j, k] * dΓ[k, 1, j] / met.J
        Tj[1, 3] += Γ[i, 1] * lc[i, j, k] * dΓ[k, 3, j] / met.J + Γ[i, 3] * lc[i, j, k] * dΓ[k, 1, j] / met.J

        

        Tj[2, 1] += Γ[i, 2] * lc[i, j, k] * dΓ[k, 1, j] / met.J + Γ[i, 1] * lc[i, j, k] * dΓ[k, 2, j] / met.J 
        #need to be multiplied by 2, flipping them reduces this to zero? unless thats what we want????
        #ie lc tensor flips.
        Tj[2, 2] += Γ[i, 2] * lc[i, j, k] * dΓ[k, 2, j] / met.J * 2 #+ Γ[k, 2] * lc[i, j, k] * dΓ[i, 2, j] / met.J
        Tj[2, 3] += Γ[i, 2] * lc[i, j, k] * dΓ[k, 3, j] / met.J + Γ[i, 3] * lc[i, j, k] * dΓ[k, 2, j] / met.J


        Tj[3, 1] += Γ[i, 3] * lc[i, j, k] * dΓ[k, 1, j] / met.J + Γ[i, 1] * lc[i, j, k] * dΓ[k, 3, j] / met.J 
        Tj[3, 2] += Γ[i, 3] * lc[i, j, k] * dΓ[k, 2, j] / met.J + Γ[i, 2] * lc[i, j, k] * dΓ[k, 3, j] / met.J
        Tj[3, 3] += Γ[i, 3] * lc[i, j, k] * dΓ[k, 3, j] / met.J * 2 #+ Γ[k, 3] * lc[i, j, k] * dΓ[i, 3, j] / met.J 

    end

    #note there should be no T[4, 4] or T[>3, >3] contribution, I think
    #looks like it is never possible to have second derive product.
    #this may make it much easier to consider every term.


    for i in 1:3, k in 1:3

        #note this can only happen once, despite there being two terms, as the second term will never have Φ_ss.
        #usure tbh, think the way we have done the double term is probably wrong
        #will need to think...
        #Ψ_s, Φ_ss can only come from first term.
        Tj[1, 4] += Γ[i, 1] * lc[i, 1, k] * Γ[k, 1] / met.J 
        #Ψ_s, Φ_sθ
        Tj[1, 5] += Γ[i, 1] * (lc[i, 1, k] * Γ[k, 2] + lc[i, 2, k] * Γ[k, 1]) / met.J
        #Ψ_s, Φ_sζ
        Tj[1, 6] += Γ[i, 1] * (lc[i, 1, k] * Γ[k, 3] + lc[i, 3, k] * Γ[k, 1]) / met.J
        #Ψ_s, Φ_θθ
        Tj[1, 7] += Γ[i, 1] * lc[i, 2, k] * Γ[k, 2] / met.J
        #Ψ_s, Φ_θζ
        Tj[1, 8] += Γ[i, 1] * (lc[i, 2, k] * Γ[k, 3] + lc[i, 3, k] * Γ[k, 2]) / met.J
        #Ψ_s, Φ_ζζ
        Tj[1, 9] += Γ[i, 1] * lc[i, 3, k] * Γ[k, 3] / met.J

        #Ψ_θ, Φ_ss
        Tj[2, 4] += Γ[i, 2] * lc[i, 1, k] * Γ[k, 1] / met.J 
        #Ψ_θ, Φ_sθ
        Tj[2, 5] += Γ[i, 2] * (lc[i, 1, k] * Γ[k, 2] + lc[i, 2, k] * Γ[k, 1]) / met.J
        #Ψ_θ, Φ_sζ
        Tj[2, 6] += Γ[i, 2] * (lc[i, 1, k] * Γ[k, 3] + lc[i, 3, k] * Γ[k, 1]) / met.J
        #Ψ_θ, Φ_θθ
        Tj[2, 7] += Γ[i, 2] * lc[i, 2, k] * Γ[k, 2] / met.J
        #Ψ_θ, Φ_θζ
        Tj[2, 8] += Γ[i, 2] * (lc[i, 2, k] * Γ[k, 3] + lc[i, 3, k] * Γ[k, 2]) / met.J
        #Ψ_θ, Φ_ζζ
        Tj[2, 9] += Γ[i, 2] * lc[i, 3, k] * Γ[k, 3] / met.J

        #Ψ_ζ, Φ_ss
        Tj[3, 4] += Γ[i, 3] * lc[i, 1, k] * Γ[k, 1] / met.J 
        #Ψ_ζ, Φ_sθ
        Tj[3, 5] += Γ[i, 3] * (lc[i, 1, k] * Γ[k, 2] + lc[i, 2, k] * Γ[k, 1]) / met.J
        #Ψ_ζ, Φ_sζ
        Tj[3, 6] += Γ[i, 3] * (lc[i, 1, k] * Γ[k, 3] + lc[i, 3, k] * Γ[k, 1]) / met.J
        #Ψ_ζ, Φ_θθ
        Tj[3, 7] += Γ[i, 3] * lc[i, 2, k] * Γ[k, 2] / met.J
        #Ψ_ζ, Φ_θζ
        Tj[3, 8] += Γ[i, 3] * (lc[i, 2, k] * Γ[k, 3] + lc[i, 3, k] * Γ[k, 2]) / met.J
        #Ψ_ζ, Φ_ζζ
        Tj[3, 9] += Γ[i, 3] * lc[i, 3, k] * Γ[k, 3] / met.J


        #now it becomes more complicated!
        #Ψ_ss, Φ_s
        Tj[4, 1] += Γ[k, 1] * lc[i, 1, k] * Γ[i, 1] / met.J
        #Ψ_ss, Φ_θ
        Tj[4, 2] += Γ[k, 1] * lc[i, 1, k] * Γ[i, 2] / met.J
        #Ψ_ss, Φ_ζ
        Tj[4, 3] += Γ[k, 1] * lc[i, 1, k] * Γ[i, 3] / met.J

        #Ψ_sθ, Φ_s
        Tj[5, 1] += (Γ[k, 1] * lc[i, 2, k] + Γ[k, 2]*lc[i, 1, k]) * Γ[i, 1] / met.J
        #Ψ_sθ, Φ_θ
        Tj[5, 2] += (Γ[k, 1] * lc[i, 2, k] + Γ[k, 2]*lc[i, 1, k]) * Γ[i, 2] / met.J
        #Ψ_sθ, Φ_ζ
        Tj[5, 3] += (Γ[k, 1] * lc[i, 2, k] + Γ[k, 2]*lc[i, 1, k]) * Γ[i, 3] / met.J

        #Ψ_sζ, Φ_s
        Tj[6, 1] += (Γ[k, 1] * lc[i, 3, k] + Γ[k, 3]*lc[i, 1, k]) * Γ[i, 1] / met.J
        #Ψ_sζ, Φ_θ
        Tj[6, 2] += (Γ[k, 1] * lc[i, 3, k] + Γ[k, 3]*lc[i, 1, k]) * Γ[i, 2] / met.J
        #Ψ_sζ, Φ_ζ
        Tj[6, 3] += (Γ[k, 1] * lc[i, 3, k] + Γ[k, 3]*lc[i, 1, k]) * Γ[i, 3] / met.J

        #Ψ_θθ, Φ_s
        Tj[7, 1] += Γ[k, 2] * lc[i, 2, k] * Γ[i, 1] / met.J
        #Ψ_θθ, Φ_θ
        Tj[7, 2] += Γ[k, 2] * lc[i, 2, k] * Γ[i, 2] / met.J
        #Ψ_θθ, Φ_ζ
        Tj[7, 3] += Γ[k, 2] * lc[i, 2, k] * Γ[i, 3] / met.J
  
        #Ψ_θζ, Φ_s
        Tj[8, 1] += (Γ[k, 2] * lc[i, 3, k] + Γ[k, 3]*lc[i, 2, k]) * Γ[i, 1] / met.J
        #Ψ_θζ, Φ_θ
        Tj[8, 2] += (Γ[k, 2] * lc[i, 3, k] + Γ[k, 3]*lc[i, 2, k]) * Γ[i, 2] / met.J
        #Ψ_θζ, Φ_ζ
        Tj[8, 3] += (Γ[k, 2] * lc[i, 3, k] + Γ[k, 3]*lc[i, 2, k]) * Γ[i, 3] / met.J

        #Ψ_ζζ, Φ_s
        Tj[9, 1] += Γ[k, 3] * lc[i, 3, k] * Γ[i, 1] / met.J
        #Ψ_ζζ, Φ_θ
        Tj[9, 2] += Γ[k, 3] * lc[i, 3, k] * Γ[i, 2] / met.J
        #Ψ_ζζ, Φ_ζ
        Tj[9, 3] += Γ[k, 3] * lc[i, 3, k] * Γ[i, 3] / met.J
    end

    #display(Tj[1:3, 1:3])
    return Tj

end

#somehow, this gives the same form of W as before, pretty wild...
#probably going to use the old form lol.
function new_compute_Tl(met, B)

    #makes more sense for this to act on W straight away I think.

    #our first step is to define the matrix C, such that
    #C @ x = ∂_i (b^j∂_j Φ)
    # = ∂_s (b^sΦ_s + b^θ Φ_θ + b^ζ Φ_ζ) + ∂_θ (b^sΦ_s + b^θ Φ_θ + b^ζ Φ_ζ) + ∂_ζ (b^sΦ_s + b^θ Φ_θ + b^ζ Φ_ζ)
    #with x = (Φ_s, Φ_θ, Φ_ζ, Φ_ss, Φ_sθ, Φ_sζ, Φ_θθ, Φ_θζ, Φ_ζζ)

    #lets try and do this very un-elegantly, just hard code each combo!
    #ie hard code what Ψ_s and Φ_s contribution will be, and repeat.

    W = zeros(Float64, 9, 9)

    W[1, 1] = 0

    for i in 1:3
        for k in 1:3
            W[1, 1] += B.db[1, i] * (met.gu[i, k] - B.b[i]*B.b[k]) * B.db[1, k]
        end
    end

    W[1, 2] = 0
    for i in 1:3
        for k in 1:3
            W[1, 2] += B.db[1, i] * (met.gu[i, k] - B.b[i]*B.b[k]) * B.db[2, k]
        end
    end

    W[1, 3] = 0
    for i in 1:3
        for k in 1:3
            W[1, 3] += B.db[1, i] * (met.gu[i, k] - B.b[i]*B.b[k]) * B.db[3, k]
        end
    end

    ################################
    W[1, 4] = 0

    for i in 1:3
        W[1, 4] += B.db[1, i] * (met.gu[i, 1] - B.b[i]*B.b[1]) * B.b[1]
    end
    
    W[1, 5] = 0

    for i in 1:3
        W[1, 5] += B.db[1, i] * (met.gu[i, 1] - B.b[i]*B.b[1]) * B.b[2] + B.db[1, i] * (met.gu[i, 2] - B.b[i]*B.b[2]) * B.b[1]
    end

    W[1, 6] = 0

    for i in 1:3
        W[1, 6] += B.db[1, i] * (met.gu[i, 1] - B.b[i]*B.b[1]) * B.b[3] + B.db[1, i] * (met.gu[i, 3] - B.b[i]*B.b[3]) * B.b[1]
    end

    W[1, 7] = 0

    for i in 1:3
        W[1, 7] += B.db[1, i] * (met.gu[i, 2] - B.b[i]*B.b[2]) * B.b[2] 
    end

    W[1, 8] = 0

    for i in 1:3
        W[1, 8] += B.db[1, i] * (met.gu[i, 2] - B.b[i]*B.b[2]) * B.b[3] + B.db[1, i] * (met.gu[i, 3] - B.b[i]*B.b[3]) * B.b[2]
    end

    W[1, 9] = 0

    for i in 1:3
        W[1, 9] += B.db[1, i] * (met.gu[i, 3] - B.b[i]*B.b[3]) * B.b[3] 
    end

    ##############

    W[2, 1] = 0

    for i in 1:3
        for k in 1:3
            W[2, 1] += B.db[2, i] * (met.gu[i, k] - B.b[i]*B.b[k]) * B.db[1, k]
        end
    end

    #Ψ_θ, Φ_θ
    W[2, 2] = 0
    for i in 1:3
        for k in 1:3
            W[2, 2] += B.db[2, i] * (met.gu[i, k] - B.b[i]*B.b[k]) * B.db[2, k]
        end
    end

    W[2, 3] = 0
    for i in 1:3
        for k in 1:3
            W[2, 3] += B.db[2, i] * (met.gu[i, k] - B.b[i]*B.b[k]) * B.db[3, k]
        end
    end

    ########################

    W[2, 4] = 0

    for i in 1:3
        W[2, 4] += B.db[2, i] * (met.gu[i, 1] - B.b[i]*B.b[1]) * B.b[1]
    end
    
    W[2, 5] = 0

    for i in 1:3
        W[2, 5] += B.db[2, i] * (met.gu[i, 1] - B.b[i]*B.b[1]) * B.b[2] + B.db[2, i] * (met.gu[i, 2] - B.b[i]*B.b[2]) * B.b[1]
    end

    W[2, 6] = 0

    for i in 1:3
        W[2, 6] += B.db[2, i] * (met.gu[i, 1] - B.b[i]*B.b[1]) * B.b[3] + B.db[2, i] * (met.gu[i, 3] - B.b[i]*B.b[3]) * B.b[1]
    end

    W[2, 7] = 0

    for i in 1:3
        W[2, 7] += B.db[2, i] * (met.gu[i, 2] - B.b[i]*B.b[2]) * B.b[2] 
    end

    W[2, 8] = 0

    for i in 1:3
        W[2, 8] += B.db[2, i] * (met.gu[i, 2] - B.b[i]*B.b[2]) * B.b[3] + B.db[2, i] * (met.gu[i, 3] - B.b[i]*B.b[3]) * B.b[2]
    end

    W[2, 9] = 0

    for i in 1:3
        W[2, 9] += B.db[2, i] * (met.gu[i, 3] - B.b[i]*B.b[3]) * B.b[3] 
    end

    #############

    W[3, 1] = 0

    for i in 1:3
        for k in 1:3
            W[3, 1] += B.db[3, i] * (met.gu[i, k] - B.b[i]*B.b[k]) * B.db[1, k]
        end
    end

    W[3, 2] = 0
    for i in 1:3
        for k in 1:3
            W[3, 2] += B.db[3, i] * (met.gu[i, k] - B.b[i]*B.b[k]) * B.db[2, k]
        end
    end

    W[3, 3] = 0
    for i in 1:3
        for k in 1:3
            W[3, 3] += B.db[3, i] * (met.gu[i, k] - B.b[i]*B.b[k]) * B.db[3, k]
        end
    end

    

    ###########################

    W[3, 4] = 0

    for i in 1:3
        W[3, 4] += B.db[3, i] * (met.gu[i, 1] - B.b[i]*B.b[1]) * B.b[1]
    end
    
    W[3, 5] = 0

    for i in 1:3
        W[3, 5] += B.db[3, i] * (met.gu[i, 1] - B.b[i]*B.b[1]) * B.b[2] + B.db[3, i] * (met.gu[i, 2] - B.b[i]*B.b[2]) * B.b[1]
    end

    W[3, 6] = 0

    for i in 1:3
        W[3, 6] += B.db[3, i] * (met.gu[i, 1] - B.b[i]*B.b[1]) * B.b[3] + B.db[3, i] * (met.gu[i, 3] - B.b[i]*B.b[3]) * B.b[1]
    end

    W[3, 7] = 0

    for i in 1:3
        W[3, 7] += B.db[3, i] * (met.gu[i, 2] - B.b[i]*B.b[2]) * B.b[2] 
    end

    W[3, 8] = 0

    for i in 1:3
        W[3, 8] += B.db[3, i] * (met.gu[i, 2] - B.b[i]*B.b[2]) * B.b[3] + B.db[3, i] * (met.gu[i, 3] - B.b[i]*B.b[3]) * B.b[2]
    end

    W[3, 9] = 0

    for i in 1:3
        W[3, 9] += B.db[3, i] * (met.gu[i, 3] - B.b[i]*B.b[3]) * B.b[3] 
    end



    #############################

    W[4, 1] = 0

    for k in 1:3
        W[4, 1] += B.b[1] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[1, k]
    end

    W[4, 2] = 0

    for k in 1:3
        W[4, 2] += B.b[1] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[2, k]
    end

    W[4, 3] = 0

    for k in 1:3
        W[4, 3] += B.b[1] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[3, k]
    end

    ############

    W[4, 4] = 0

    #for k in 1:3
    W[4, 4] += B.b[1] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[1]
    #end

    W[4, 5] = 0

    #for k in 1:3
    W[4, 5] += B.b[1] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[2] + B.b[1] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[1]
    #end

    W[4, 6] = 0

    #for k in 1:3
    W[4, 6] += B.b[1] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[3] + B.b[1] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[1]


    W[4, 7] = 0

    #for k in 1:3
    W[4, 7] += B.b[1] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[2] 

    W[4, 8] = 0

    #for k in 1:3
    W[4, 8] += B.b[1] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[3] + B.b[1] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[2]

    W[4, 9] = 0

    #for k in 1:3
    W[4, 9] += B.b[1] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[3] 


    #####################

   W[5, 1] = 0

    for k in 1:3
        W[5, 1] += B.b[1] * (met.gu[2, k] - B.b[2]*B.b[k]) * B.db[1, k] + B.b[2] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[1, k]
    end

    W[5, 2] = 0

    for k in 1:3
        W[5, 2] += B.b[1] * (met.gu[2, k] - B.b[2]*B.b[k]) * B.db[2, k] + B.b[2] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[2, k]
    end

    W[5, 3] = 0

    for k in 1:3
        W[5, 3] += B.b[1] * (met.gu[2, k] - B.b[2]*B.b[k]) * B.db[3, k] + B.b[2] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[3, k]
    end

    ###############

    W[5, 4] = 0

    #for k in 1:3
    W[5, 4] += B.b[1] * (met.gu[2, 1] - B.b[2]*B.b[1]) * B.b[1] + B.b[2] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[1]
    #end

    W[5, 5] = 0

    #for k in 1:3
    W[5, 5] += (B.b[1] * (met.gu[2, 1] - B.b[2]*B.b[1]) * B.b[2] + B.b[1] * (met.gu[2, 2] - B.b[2]*B.b[2]) * B.b[1] 
                + B.b[2] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[2] + B.b[2] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[1])
    #end

    W[5, 6] = 0

    W[5, 6] += (B.b[1] * (met.gu[2, 1] - B.b[2]*B.b[1]) * B.b[3] + B.b[1] * (met.gu[2, 3] - B.b[2]*B.b[3]) * B.b[1] 
                + B.b[2] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[3] + B.b[2] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[1])
    

    W[5, 7] = 0

    #for k in 1:3
    W[5, 7] += B.b[1] * (met.gu[2, 2] - B.b[2]*B.b[2]) * B.b[2] + B.b[2] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[2] 

    W[5, 8] = 0

    #for k in 1:3
    W[5, 8] += (B.b[1] * (met.gu[2, 2] - B.b[2]*B.b[2]) * B.b[3] + B.b[1] * (met.gu[2, 3] - B.b[2]*B.b[3]) * B.b[2]
                + B.b[2] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[3] + B.b[2] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[2])

    W[5, 9] = 0

    #for k in 1:3
    W[5, 9] += B.b[1] * (met.gu[2, 3] - B.b[2]*B.b[3]) * B.b[3] + B.b[2] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[3] 


    #######################

    W[6, 1] = 0

    for k in 1:3
        W[6, 1] += B.b[1] * (met.gu[3, k] - B.b[3]*B.b[k]) * B.db[1, k] + B.b[3] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[1, k]
    end

    W[6, 2] = 0

    for k in 1:3
        W[6, 2] += B.b[1] * (met.gu[3, k] - B.b[3]*B.b[k]) * B.db[2, k] + B.b[3] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[2, k]
    end

    W[6, 3] = 0

    for k in 1:3
        W[6, 3] += B.b[1] * (met.gu[3, k] - B.b[3]*B.b[k]) * B.db[3, k] + B.b[3] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[3, k]
    end

    ###############

    W[6, 4] = 0

    #for k in 1:3
    W[6, 4] += B.b[1] * (met.gu[3, 1] - B.b[3]*B.b[1]) * B.b[1] + B.b[3] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[1]
    #end

    W[6, 5] = 0

    #for k in 1:3
    W[6, 5] += (B.b[1] * (met.gu[3, 1] - B.b[3]*B.b[1]) * B.b[2] + B.b[1] * (met.gu[3, 2] - B.b[3]*B.b[2]) * B.b[1] 
                + B.b[3] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[2] + B.b[3] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[1])
    #end

    W[6, 6] = 0

    W[6, 6] += (B.b[1] * (met.gu[3, 1] - B.b[3]*B.b[1]) * B.b[3] + B.b[1] * (met.gu[3, 3] - B.b[3]*B.b[3]) * B.b[1] 
                + B.b[3] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[3] + B.b[3] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[1])
    

    W[6, 7] = 0

    #for k in 1:3
    W[6, 7] += B.b[1] * (met.gu[3, 2] - B.b[3]*B.b[2]) * B.b[2] + B.b[3] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[2] 

    W[6, 8] = 0

    #for k in 1:3
    W[6, 8] += (B.b[1] * (met.gu[3, 2] - B.b[3]*B.b[2]) * B.b[3] + B.b[1] * (met.gu[3, 3] - B.b[3]*B.b[3]) * B.b[2]
                + B.b[3] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[3] + B.b[3] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[2])

    W[6, 9] = 0

    #for k in 1:3
    W[6, 9] += B.b[1] * (met.gu[3, 3] - B.b[3]*B.b[3]) * B.b[3] + B.b[3] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[3] 

    ######################

    W[7, 1] = 0

    for k in 1:3
        W[7, 1] += B.b[2] * (met.gu[2, k] - B.b[2]*B.b[k]) * B.db[1, k] #+ B.b[2] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[1, k]
    end

    W[7, 2] = 0

    for k in 1:3
        W[7, 2] += B.b[2] * (met.gu[2, k] - B.b[2]*B.b[k]) * B.db[2, k] #+ B.b[2] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[2, k]
    end

    W[7, 3] = 0

    for k in 1:3
        W[7, 3] += B.b[2] * (met.gu[2, k] - B.b[2]*B.b[k]) * B.db[3, k] #+ B.b[2] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[3, k]
    end

    ###############

    W[7, 4] = 0

    #for k in 1:3
    W[7, 4] += B.b[2] * (met.gu[2, 1] - B.b[2]*B.b[1]) * B.b[1] #+ B.b[2] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[1]
    #end

    W[7, 5] = 0

    #for k in 1:3
    W[7, 5] += (B.b[2] * (met.gu[2, 1] - B.b[2]*B.b[1]) * B.b[2] + B.b[2] * (met.gu[2, 2] - B.b[2]*B.b[2]) * B.b[1] )
                #+ B.b[2] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[2] + B.b[2] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[1])
    #end

    W[7, 6] = 0

    W[7, 6] += (B.b[2] * (met.gu[2, 1] - B.b[2]*B.b[1]) * B.b[3] + B.b[2] * (met.gu[2, 3] - B.b[2]*B.b[3]) * B.b[1]) 
                #+ B.b[2] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[3] + B.b[2] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[1])
    

    W[7, 7] = 0

    #for k in 1:3
    W[7, 7] += B.b[2] * (met.gu[2, 2] - B.b[2]*B.b[2]) * B.b[2] #+ B.b[2] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[2] 

    W[7, 8] = 0

    #for k in 1:3
    W[7, 8] += (B.b[2] * (met.gu[2, 2] - B.b[2]*B.b[2]) * B.b[3] + B.b[2] * (met.gu[2, 3] - B.b[2]*B.b[3]) * B.b[2])
                #+ B.b[2] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[3] + B.b[2] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[2])

    W[7, 9] = 0

    #for k in 1:3
    W[7, 9] += B.b[2] * (met.gu[2, 3] - B.b[2]*B.b[3]) * B.b[3] #+ B.b[2] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[3] 


    ##########################

    W[8, 1] = 0

    for k in 1:3
        W[8, 1] += B.b[2] * (met.gu[3, k] - B.b[3]*B.b[k]) * B.db[1, k] + B.b[3] * (met.gu[2, k] - B.b[2]*B.b[k]) * B.db[1, k]
    end

    W[8, 2] = 0

    for k in 1:3
        W[8, 2] += B.b[2] * (met.gu[3, k] - B.b[3]*B.b[k]) * B.db[2, k] + B.b[3] * (met.gu[2, k] - B.b[2]*B.b[k]) * B.db[2, k]
    end

    W[8, 3] = 0

    for k in 1:3
        W[8, 3] += B.b[2] * (met.gu[3, k] - B.b[3]*B.b[k]) * B.db[3, k] + B.b[3] * (met.gu[2, k] - B.b[2]*B.b[k]) * B.db[3, k]
    end

    ###############

    W[8, 4] = 0

    #for k in 1:3
    W[8, 4] += B.b[2] * (met.gu[3, 1] - B.b[3]*B.b[1]) * B.b[1] + B.b[3] * (met.gu[2, 1] - B.b[2]*B.b[1]) * B.b[1]
    #end

    W[8, 5] = 0

    #for k in 1:3
    W[8, 5] += (B.b[2] * (met.gu[3, 1] - B.b[3]*B.b[1]) * B.b[2] + B.b[2] * (met.gu[3, 2] - B.b[3]*B.b[2]) * B.b[1] 
                + B.b[3] * (met.gu[2, 1] - B.b[2]*B.b[1]) * B.b[2] + B.b[3] * (met.gu[2, 2] - B.b[2]*B.b[2]) * B.b[1])
    #end

    W[8, 6] = 0

    W[8, 6] += (B.b[2] * (met.gu[3, 1] - B.b[3]*B.b[1]) * B.b[3] + B.b[2] * (met.gu[3, 3] - B.b[3]*B.b[3]) * B.b[1] 
                + B.b[3] * (met.gu[2, 1] - B.b[2]*B.b[1]) * B.b[3] + B.b[3] * (met.gu[2, 3] - B.b[2]*B.b[3]) * B.b[1])
    

    W[8, 7] = 0

    #for k in 1:3
    W[8, 7] += B.b[2] * (met.gu[3, 2] - B.b[3]*B.b[2]) * B.b[2] + B.b[3] * (met.gu[2, 2] - B.b[2]*B.b[2]) * B.b[2] 

    W[8, 8] = 0

    #for k in 1:3
    W[8, 8] += (B.b[2] * (met.gu[3, 2] - B.b[3]*B.b[2]) * B.b[3] + B.b[2] * (met.gu[3, 3] - B.b[3]*B.b[3]) * B.b[2]
                + B.b[3] * (met.gu[2, 2] - B.b[2]*B.b[2]) * B.b[3] + B.b[3] * (met.gu[2, 3] - B.b[2]*B.b[3]) * B.b[2])

    W[8, 9] = 0

    #for k in 1:3
    W[8, 9] += B.b[2] * (met.gu[3, 3] - B.b[3]*B.b[3]) * B.b[3] + B.b[3] * (met.gu[2, 3] - B.b[2]*B.b[3]) * B.b[3] 


    #######################################
    W[9, 1] = 0

    for k in 1:3
        W[9, 1] += B.b[3] * (met.gu[3, k] - B.b[3]*B.b[k]) * B.db[1, k] #+ B.b[2] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[1, k]
    end

    W[9, 2] = 0

    for k in 1:3
        W[9, 2] += B.b[3] * (met.gu[3, k] - B.b[3]*B.b[k]) * B.db[2, k] #+ B.b[2] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[2, k]
    end

    W[9, 3] = 0

    for k in 1:3
        W[9, 3] += B.b[3] * (met.gu[3, k] - B.b[3]*B.b[k]) * B.db[3, k] #+ B.b[2] * (met.gu[1, k] - B.b[1]*B.b[k]) * B.db[3, k]
    end

    ###############

    W[9, 4] = 0

    #for k in 1:3
    W[9, 4] += B.b[3] * (met.gu[3, 1] - B.b[3]*B.b[1]) * B.b[1] #+ B.b[2] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[1]
    #end

    W[9, 5] = 0

    #for k in 1:3
    W[9, 5] += (B.b[3] * (met.gu[3, 1] - B.b[3]*B.b[1]) * B.b[2] + B.b[3] * (met.gu[3, 2] - B.b[3]*B.b[2]) * B.b[1] )
                #+ B.b[2] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[2] + B.b[2] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[1])
    #end

    W[9, 6] = 0

    W[9, 6] += (B.b[3] * (met.gu[3, 1] - B.b[3]*B.b[1]) * B.b[3] + B.b[3] * (met.gu[3, 3] - B.b[3]*B.b[3]) * B.b[1] )
                #+ B.b[2] * (met.gu[1, 1] - B.b[1]*B.b[1]) * B.b[3] + B.b[2] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[1])
    

    W[9, 7] = 0

    #for k in 1:3
    W[9, 7] += B.b[3] * (met.gu[3, 2] - B.b[3]*B.b[2]) * B.b[2] #+ B.b[2] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[2] 

    W[9, 8] = 0

    #for k in 1:3
    W[9, 8] += (B.b[3] * (met.gu[3, 2] - B.b[3]*B.b[2]) * B.b[3] + B.b[3] * (met.gu[3, 3] - B.b[3]*B.b[3]) * B.b[2])
                #+ B.b[2] * (met.gu[1, 2] - B.b[1]*B.b[2]) * B.b[3] + B.b[2] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[2])

    W[9, 9] = 0

    #for k in 1:3
    W[9, 9] += B.b[3] * (met.gu[3, 3] - B.b[3]*B.b[3]) * B.b[3] #+ B.b[2] * (met.gu[1, 3] - B.b[1]*B.b[3]) * B.b[3] 


    return W

end