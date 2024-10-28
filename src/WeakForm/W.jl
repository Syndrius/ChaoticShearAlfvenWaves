
#this is the worst file in the whole package!!

"""
    compute_W!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, met::MetT, B::BFieldT)

Computes the W matrix for the weak form at a single coordinate. Awful function that is built from awful functions! This file needs work!
"""
function compute_W!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, met::MetT, B::BFieldT, n::Float64, ωcap2::Float64)

    #now we want to combine both Tj and Tl into one, this will be cooked!
    #Tl = zeros(9, 9)
    #Tj = zeros(3, 9)
    #display("og")

    #contributions to W are split into two.
    W[:, :] = compute_Tl(met, B) .* met.J #Tl is fine I think.

    #for islands we assume current is negligible.
    W[:, :] -= compute_Tj(met, B) .* met.J .* jparonB(met, B) ./ 2


    #W[:, :] -= compute_cap(met, B)

    

    for j=1:3, i=1:3
        W[i, j] += ωcap2 * n*(met.gu[i, j] - B.b[i]*B.b[j]) * met.J / B.mag_B^2
    end
    #Tj = compute_Tj(met, B)

    #display(W[:, :, co...])
    #need to set equal here, the fill function is not working as intended.
    

    #this could be total garbage.
    
    #for μ=1:9, i=1:3
    #    W[i, μ] += Tj[i, μ] * met.J
    #    W[μ, i] += Tj[i, μ] * met.J
    #end 
    
    
    #println(W)

    
end


#think this is actually not needed
#as island case does not compute dB or dg, so current term will always be zero regardless.
#also, we may eventually want those additions,
#however, this should be a bit faster.
function compute_isl_W!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, met::MetT, B::BFieldT, n::Float64)

    #now we want to combine both Tj and Tl into one, this will be cooked!
    #Tl = zeros(9, 9)
    #Tj = zeros(3, 9)
    #display("og")

    #contributions to W are split into two.

    W[:, :] = compute_Tl(met, B) .* met.J #Tl is fine I think.

    
    
end


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

#const lct = cat(3, [0 0 0; 0 0 1; 0 -1 0], [0 0 -1; 0 0 0; 1 0 0], [0 1 0; -1 0 0; 0 0 0], dims=(3, 3, 3))
"""
Function for the levi-civita tensor.
"""
function get_lc_tensor()
    #probably a more elegant way to do this.
    lc = zeros(3, 3, 3)
    lc[1, 2, 3] = 1
    lc[2, 3, 1] = 1
    lc[3, 1, 2] = 1
    lc[3, 2, 1] = -1
    lc[1, 3, 2] = -1
    lc[2, 1, 3] = -1
    return lc
end


const lct = get_lc_tensor()

"""
    jparonB(met::MetT, B::BFieldT)

Computes the parrallel current divided by the magnitude of B, needed for the current term of W.
"""
function jparonB(met::MetT, B::BFieldT)

    J = zeros(3)

    #J^i = (∇×B)^i = 1/J * ϵ^{ijk}∂_j B_k = 1/J * ϵ^{ijk}∂_j (g_{kl} B^l)

    for i in 1:3, j in 1:3, k in 1:3
        J[i] += 1 / met.J * lct[i, j, k] * (dot(met.gl[k, :], B.dB[:, j]) + dot(met.dgl[k, :, j], B.B[:]))
    end
    jpar = 0
    for i in 1:3, j in 1:3
        jpar += met.gl[i, j] * B.b[i] * J[j]
    end
    return jpar/B.mag_B
end

function dδ(i::Int64, j::Int64)
    if i==j
        return 1.0
    else
        return 0.0
    end
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