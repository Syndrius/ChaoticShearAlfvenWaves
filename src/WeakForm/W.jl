
#this is the worst file in the whole package!!

"""
Computes the W matrix for the weak form at a single coordinate. Awful function that is built from awful functions! This file needs work!

# Args 
- I Stupid type due to using views.
- met::MetT Struct containing the metric information.
- B::BfieldT Struct containing the magnetic field information.
"""
function compute_W!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, met::MetT, B::BFieldT)

    #now we want to combine both Tj and Tl into one, this will be cooked!
    #Tl = zeros(9, 9)
    #Tj = zeros(3, 9)
    #display("og")

    #contributions to W are split into two.
    W[:, :] = compute_Tl(met, B) .* met.J
    
    Tj = compute_Tj(met, B)

    #display(W[:, :, co...])
    #need to set equal here, the fill function is not working as intended.
    

    #this could be total garbage.
    
    for μ=1:9, i=1:3
        W[i, μ] += Tj[i, μ] * met.J
        W[μ, i] += Tj[i, μ] * met.J
    end 
    
    
    #println(W)

    
end



#this might be the worst function ever written...
function compute_Tj(met::MetT, B::BFieldT)

    Γ = zeros(3, 3)
    dΓ = zeros(3, 3, 3)

    Tj = zeros(3, 9)

    """
    for i in 1:3
        for j in 1:3
            for m in 1:3
                #Γ^j_i = (δ_i^j - g_{im} b^m b^j)
                #fk me we are adding the dirac term 3 times...
                #how was this ever working????
                Γ[i, j] += dirac_delta(i, j) - met.gl[i, m]*B.b[m]*B.b[j]

                for k in 1:3
                    #∂_k (Γ^j_i) = - b^m b^j ∂_k(g_{im}) - g_{im} b^m ∂_k(b^j) - g_{im} b^j ∂_k(b^m)
                    dΓ[i, j, k] -= (B.b[m] * B.b[j] * met.dgl[i, m, k] +
                                    B.b[m] * met.gl[i, m] * B.db[j, k] + 
                                    B.b[j] * met.gl[i, m] * B.db[m, k])
                end
            end
        end
    end

    """


    #first we deal with the diract delta!
    #this is julia identity, may need a clearer solution given we use an I matrix elsewhere!
    #this is different to the other case!!!! wtf! Nope I am just a goose
    #display(I)
    #Γ += I
    
    #poor solution for adding the dirac delta contribution
    Γ[1, 1] = 1
    Γ[2, 2] = 1
    Γ[3, 3] = 1
    #display(Γ)


    for i in 1:3
        #this is better but a fkn garbage way to do this!
        Γ[1, 1] -= met.gl[1, i] * B.b[i] * B.b[1]
        Γ[1, 2] -= met.gl[1, i] * B.b[i] * B.b[2]
        Γ[1, 3] -= met.gl[1, i] * B.b[i] * B.b[3]

        Γ[2, 1] -= met.gl[2, i] * B.b[i] * B.b[1]
        Γ[2, 2] -= met.gl[2, i] * B.b[i] * B.b[2]
        Γ[2, 3] -= met.gl[2, i] * B.b[i] * B.b[3]

        Γ[3, 1] -= met.gl[3, i] * B.b[i] * B.b[1]
        Γ[3, 2] -= met.gl[3, i] * B.b[i] * B.b[2]
        Γ[3, 3] -= met.gl[3, i] * B.b[i] * B.b[3]

        #awful way of doing this but should be more reliable
        dΓ[1, 1, 1] -= (met.dgl[1, i, 1] * B.b[i] * B.b[1] 
                    + met.gl[1, i] * B.db[i, 1] * B.b[1]
                    + met.gl[1, i] * B.b[i] * B.db[1, 1])
        dΓ[1, 1, 2] -= (met.dgl[1, i, 2] * B.b[i] * B.b[1] 
                    + met.gl[1, i] * B.db[i, 2] * B.b[1]
                    + met.gl[1, i] * B.b[i] * B.db[1, 2])
        dΓ[1, 1, 3] -= (met.dgl[1, i, 3] * B.b[i] * B.b[1] 
                    + met.gl[1, i] * B.db[i, 3] * B.b[1]
                    + met.gl[1, i] * B.b[i] * B.db[1, 3])
        

        dΓ[1, 2, 1] -= (met.dgl[1, i, 1] * B.b[i] * B.b[2] 
                    + met.gl[1, i] * B.db[i, 1] * B.b[2]
                    + met.gl[1, i] * B.b[i] * B.db[2, 1])
        dΓ[1, 2, 2] -= (met.dgl[1, i, 2] * B.b[i] * B.b[2] 
                    + met.gl[1, i] * B.db[i, 2] * B.b[2]
                    + met.gl[1, i] * B.b[i] * B.db[2, 2])
        dΓ[1, 2, 3] -= (met.dgl[1, i, 3] * B.b[i] * B.b[2] 
                    + met.gl[1, i] * B.db[i, 3] * B.b[2]
                    + met.gl[1, i] * B.b[i] * B.db[2, 3])

        dΓ[1, 3, 1] -= (met.dgl[1, i, 1] * B.b[i] * B.b[3] 
                    + met.gl[1, i] * B.db[i, 1] * B.b[3]
                    + met.gl[1, i] * B.b[i] * B.db[3, 1])
        dΓ[1, 3, 2] -= (met.dgl[1, i, 2] * B.b[i] * B.b[3] 
                    + met.gl[1, i] * B.db[i, 2] * B.b[3]
                    + met.gl[1, i] * B.b[i] * B.db[3, 2])
        dΓ[1, 3, 3] -= (met.dgl[1, i, 3] * B.b[i] * B.b[3] 
                    + met.gl[1, i] * B.db[i, 3] * B.b[3]
                    + met.gl[1, i] * B.b[i] * B.db[3, 3])
        


        dΓ[2, 1, 1] -= (met.dgl[2, i, 1] * B.b[i] * B.b[1] 
                    + met.gl[2, i] * B.db[i, 1] * B.b[1]
                    + met.gl[2, i] * B.b[i] * B.db[1, 1])
        dΓ[2, 1, 2] -= (met.dgl[2, i, 2] * B.b[i] * B.b[1] 
                    + met.gl[2, i] * B.db[i, 2] * B.b[1]
                    + met.gl[2, i] * B.b[i] * B.db[1, 2])
        dΓ[2, 1, 3] -= (met.dgl[2, i, 3] * B.b[i] * B.b[1] 
                    + met.gl[2, i] * B.db[i, 3] * B.b[1]
                    + met.gl[2, i] * B.b[i] * B.db[1, 3])
        

        dΓ[2, 2, 1] -= (met.dgl[2, i, 1] * B.b[i] * B.b[2] 
                    + met.gl[2, i] * B.db[i, 1] * B.b[2]
                    + met.gl[2, i] * B.b[i] * B.db[2, 1])
        dΓ[2, 2, 2] -= (met.dgl[2, i, 2] * B.b[i] * B.b[2] 
                    + met.gl[2, i] * B.db[i, 2] * B.b[2]
                    + met.gl[2, i] * B.b[i] * B.db[2, 2])
        dΓ[2, 2, 3] -= (met.dgl[2, i, 3] * B.b[i] * B.b[2] 
                    + met.gl[2, i] * B.db[i, 3] * B.b[2]
                    + met.gl[2, i] * B.b[i] * B.db[2, 3])

        dΓ[2, 3, 1] -= (met.dgl[2, i, 1] * B.b[i] * B.b[3] 
                    + met.gl[2, i] * B.db[i, 1] * B.b[3]
                    + met.gl[2, i] * B.b[i] * B.db[3, 1])
        dΓ[2, 3, 2] -= (met.dgl[2, i, 2] * B.b[i] * B.b[3] 
                    + met.gl[2, i] * B.db[i, 2] * B.b[3]
                    + met.gl[2, i] * B.b[i] * B.db[3, 2])
        dΓ[2, 3, 3] -= (met.dgl[2, i, 3] * B.b[i] * B.b[3] 
                    + met.gl[2, i] * B.db[i, 3] * B.b[3]
                    + met.gl[2, i] * B.b[i] * B.db[3, 3])




        dΓ[3, 1, 1] -= (met.dgl[3, i, 1] * B.b[i] * B.b[1] 
                    + met.gl[3, i] * B.db[i, 1] * B.b[1]
                    + met.gl[3, i] * B.b[i] * B.db[1, 1])
        dΓ[3, 1, 2] -= (met.dgl[3, i, 2] * B.b[i] * B.b[1] 
                    + met.gl[3, i] * B.db[i, 2] * B.b[1]
                    + met.gl[3, i] * B.b[i] * B.db[1, 2])
        dΓ[3, 1, 3] -= (met.dgl[3, i, 3] * B.b[i] * B.b[1] 
                    + met.gl[3, i] * B.db[i, 3] * B.b[1]
                    + met.gl[3, i] * B.b[i] * B.db[1, 3])
        

        dΓ[3, 2, 1] -= (met.dgl[3, i, 1] * B.b[i] * B.b[2] 
                    + met.gl[3, i] * B.db[i, 1] * B.b[2]
                    + met.gl[3, i] * B.b[i] * B.db[2, 1])
        dΓ[3, 2, 2] -= (met.dgl[3, i, 2] * B.b[i] * B.b[2] 
                    + met.gl[3, i] * B.db[i, 2] * B.b[2]
                    + met.gl[3, i] * B.b[i] * B.db[2, 2])
        dΓ[3, 2, 3] -= (met.dgl[3, i, 3] * B.b[i] * B.b[2] 
                    + met.gl[3, i] * B.db[i, 3] * B.b[2]
                    + met.gl[3, i] * B.b[i] * B.db[2, 3])

        dΓ[3, 3, 1] -= (met.dgl[3, i, 1] * B.b[i] * B.b[3] 
                    + met.gl[3, i] * B.db[i, 1] * B.b[3]
                    + met.gl[3, i] * B.b[i] * B.db[3, 1])
        dΓ[3, 3, 2] -= (met.dgl[3, i, 2] * B.b[i] * B.b[3] 
                    + met.gl[3, i] * B.db[i, 2] * B.b[3]
                    + met.gl[3, i] * B.b[i] * B.db[3, 2])
        dΓ[3, 3, 3] -= (met.dgl[3, i, 3] * B.b[i] * B.b[3] 
                    + met.gl[3, i] * B.db[i, 3] * B.b[3]
                    + met.gl[3, i] * B.b[i] * B.db[3, 3])
    end

    #display(B.b[:])
    #display(Γ)
    #display(dΓ)
    #display(dΓ[:, :, 2])
    #display(sqrt(-1))
    
    """
    for i in 1:3
        for j in 1:3
            for k in 1:3
                for n in 1:3
                    for q in 1:3

                        #think we may need to be real thorough and check all of this.
                        #first deriv part
                        Tj[n, q] += Γ[i, n] / met.J * lc[i, j, k] * dΓ[k, q, j]
                        #second deriv part
                        #this gives worse results
                        #something is either wrong with this loop or with the d_to_ind func.
                        #this is defs fked no matter what lol.
                        #Tj[n, d_to_ind(j, q)] += Γ[i, n] / met.J * lc[i, j, k] * Γ[k, q]
                    end
                end
            end
        end
    end
    """

    #this seems to have made no difference
    
    for n in 1:3

        #r q=1
        #i=1 + i=2 + i=3
        #ϵ^123 + ϵ^132
        Tj[n, 1] = Γ[1, n] * (dΓ[3, 1, 2] - dΓ[2, 1, 3]) / met.J
        #ϵ^231 + ϵ^213
        Tj[n, 1] += Γ[2, n] * (dΓ[1, 1, 3] - dΓ[3, 1, 1]) / met.J
        #ϵ^312 + ϵ^321
        Tj[n, 1] += Γ[3, n] * (dΓ[2, 1, 1] - dΓ[1, 1, 2]) / met.J

        #θ q=2
        #i=1 + i=2 + i=3
        #ϵ^123 + ϵ^132
        Tj[n, 2] = Γ[1, n] * (dΓ[3, 2, 2] - dΓ[2, 2, 3]) / met.J
        #ϵ^231 + ϵ^213
        Tj[n, 2] += Γ[2, n] * (dΓ[1, 2, 3] - dΓ[3, 2, 1]) / met.J
        #ϵ^312 + ϵ^321
        Tj[n, 2] += Γ[3, n] * (dΓ[2, 2, 1] - dΓ[1, 2, 2]) / met.J

        #ζ q=3
        #i=1 + i=2 + i=3
        #ϵ^123 + ϵ^132
        Tj[n, 3] = Γ[1, n] * (dΓ[3, 3, 2] - dΓ[2, 3, 3]) / met.J
        #ϵ^231 + ϵ^213
        Tj[n, 3] += Γ[2, n] * (dΓ[1, 3, 3] - dΓ[3, 3, 1]) / met.J
        #ϵ^312 + ϵ^321
        Tj[n, 3] += Γ[3, n] * (dΓ[2, 3, 1] - dΓ[1, 3, 2]) / met.J

    end
    

    #sqrt(-1)
    #display(Γ[:, :])

    #display(met.dgl[:, :, 1])
    
    #this is much closer but results are still not perf, not sure if this has a mistake in it
    #or Tj[1:3] also has issues, ie the above loop could just be cooked somehow.
    for n in 1:3
        #following weak form indicies!
        #rr
        #j=1, q=1
        #ie ϵ^312 + ϵ^213, i=3, k=2 + i=2, k=3
        Tj[n, 4] = (Γ[3, n] * Γ[2, 1] - Γ[2, n] * Γ[3, 1]) / met.J

        #rθ
        #first j=2, q=1
        #ϵ^123 + ϵ^321, i=1, k=3 + i=3, k=1
        Tj[n, 5] = (Γ[1, n] * Γ[3, 1] - Γ[3, n] * Γ[1, 1]) / met.J
        #then j=1, q=2
        #ϵ^312 + ϵ^213, i=3, k=2 + i=2, k=3
        Tj[n, 5] += (Γ[3, n] * Γ[2, 2] - Γ[2, n] * Γ[3, 2]) / met.J

        #n=2 case
        #Tj[2, 5] = (Γ[1, 2] * Γ[3, 1] - Γ[3, n] * Γ[1, 1]) / met.J
        #    + (Γ[3, 2] * Γ[2, 2] - Γ[2, 2] * Γ[3, 2]) / met.J

        #rζ
        #first j=3, q=1
        #ϵ^231 + ϵ^132, i=2, k=1 + i=1, k=2
        Tj[n, 6] = (Γ[2, n] * Γ[1, 1] - Γ[1, n] * Γ[2, 1]) / met.J
        #then j=1, q=3
        #ϵ^312 + ϵ^213, i=3, k=2 + i=2, k=3
        Tj[n, 6] += (Γ[3, n] * Γ[2, 3] - Γ[2, n] * Γ[3, 3]) / met.J

        #θθ
        #j=2, q=2
        #ϵ^123 + ϵ^321, i=1, k=3 + i=3, k=1
        Tj[n, 7] = (Γ[1, n] * Γ[3, 2] - Γ[3, n] * Γ[1, 2]) / met.J

        #θζ
        #first j=3, q=2
        #ϵ^231 + ϵ^132, i=2, k=1 + i=1, k=2
        Tj[n, 8] = (Γ[2, n] * Γ[1, 2] - Γ[1, n] * Γ[2, 2]) / met.J
        #then j=2, q=3
        #ϵ^123 + ϵ^321, i=1, k=3 + i=3, k=1
        Tj[n, 8] += (Γ[1, n] * Γ[3, 3] - Γ[3, n] * Γ[1, 3]) / met.J

        #ζζ
        #j=3, q=3
        #ϵ^231 + ϵ^132, i=2, k=1 + i=1, k=2
        Tj[n, 9] = (Γ[2, n] * Γ[1, 3] - Γ[1, n] * Γ[2, 3]) / met.J
    end
    
    
    #looks like the self adjoint stuff is maybe more important than we give credit for.
    #If we return this, only contract `half` of Tj, but solve with Hermitian flag we get nonsense.
    #however solving without Hermitian flag gives us identical answers.
    #this is probably the problem we had earlier.
    #even if it accepts Hermitian, and the operator is technically Hermitian, the matrix itself is not!
    #so we get garbage.!
    #return -1 .* jparonB(B, met* Tj 
    #so jparonB changes massivly with the different metric indicies.
    #display(jparonB(B, met, co))

    #only the first 3x3 changes for different mets.
    #display(Tj )
    #display(dΓ[:, :, 1])

    #display(sqrt(-1))
    #display(jparonB(B, met, co))

    #display(2 .*Tj[1:3, 1:3])
    #return Tj #for testing!
    return -1/2 .* jparonB(met, B) .* Tj 
end



function jparonB(met::MetT, B::BFieldT)

    J = zeros(3)

    #J^i = (∇×B)^i = 1/J * ϵ^{ijk}∂_j B_k = 1/J * ϵ^{ijk}∂_j (g_{kl} B^l)

    #this is an awful way of doing this...
    #this is giving different, yet better results...
    #maybe I cooked the lc tensor somehow???
    #cannot believe the other form doesn't work...
    #can't be a coincidence that terms involving lc seem to be broken.
    #this seems to not be true now, as we fixed the metric...

    #J^1 = 1/J * ϵ^{1jk}∂_j (g_{kl} B^l) = 1/J * ϵ^{123}∂_2 (g_{l3} B^l) + 1/J * ϵ^{132}∂_3 (g_{l2} B^l)
    J[1] = 1/met.J * (B.dB[1, 2] * met.gl[1, 3] + B.B[1] * met.dgl[1, 3, 2] +
    B.dB[2, 2] * met.gl[2, 3] + B.B[2] * met.dgl[2, 3, 2] + 
    B.dB[3, 2] * met.gl[3, 3] + B.B[3] * met.dgl[3, 3, 2]) 

    J[1] -= 1/met.J * (B.dB[1, 3] * met.gl[1, 2] + B.B[1] * met.dgl[1, 2, 3] +
    B.dB[2, 3] * met.gl[2, 2] + B.B[2] * met.dgl[2, 2, 3] + 
    B.dB[3, 3] * met.gl[3, 2] + B.B[3] * met.dgl[3, 2, 3]) 


    J[2] = 1/met.J * (B.dB[1, 3] * met.gl[1, 1] + B.B[1] * met.dgl[1, 1, 3] +
    B.dB[2, 3] * met.gl[2, 1] + B.B[2] * met.dgl[2, 1, 3] + 
    B.dB[3, 3] * met.gl[3, 1] + B.B[3] * met.dgl[3, 1, 3]) 

    J[2] -= 1/met.J * (B.dB[1, 1] * met.gl[1, 3] + B.B[1] * met.dgl[1, 3, 1] +
    B.dB[2, 1] * met.gl[2, 3] + B.B[2] * met.dgl[2, 3, 1] + 
    B.dB[3, 1] * met.gl[3, 3] + B.B[3] * met.dgl[3, 3, 1]) 

    J[3] = 1/met.J * (B.dB[1, 1] * met.gl[1, 2] + B.B[1] * met.dgl[1, 2, 1] +
    B.dB[2, 1] * met.gl[2, 2] + B.B[2] * met.dgl[2, 2, 1] + 
    B.dB[3, 1] * met.gl[3, 2] + B.B[3] * met.dgl[3, 2, 1]) 

    J[3] -= 1/met.J * (B.dB[1, 2] * met.gl[1, 1] + B.B[1] * met.dgl[1, 1, 2] +
    B.dB[2, 2] * met.gl[2, 1] + B.B[2] * met.dgl[2, 1, 2] + 
    B.dB[3, 2] * met.gl[3, 1] + B.B[3] * met.dgl[3, 1, 2]) 

    #with the old metric, the current is completly cooked and wrong!
    #and yet...
    #display("current")
    #display(J)

"""
    for i in 1:3
        for j in 1:3

            for k in 1:3
                for l in 1:3


                    J[i] += 1/met.J * lc[i, j, k] * (B.dB[l, j]*met.gl[k, l] + B.B[l] * met.dgl[k, l, j])
                end
            end
        end
    end
    
    """
    jparonB = 0

    for i in 1:3
        for j in 1:3
            jparonB += B.b[i] * J[j]  / B.mag_B * met.gl[i, j]
        end
    end

    #display(jparonB)
    return jparonB
end


function compute_Tl(met::MetT, B::BFieldT)

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