"""
    function Tj!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, B::BFieldT, met::MetT, Γ::Array{Float64, 2}, dΓ::Array{Float64, 3}, K::Array{Float64})

Computes the current term contribution for W.
Tj  is given by the expression;
Γ_i^n ∂_nΨ 1/J ϵ^{ijk}(Γ_k^q∂_j∂_qΦ + ∂_j(Γ_k^q)∂_qΦ) + Γ_i^n ∂_nΦ 1/J ϵ^{ijk}(Γ_k^q∂_j∂_qΨ + ∂_j(Γ_k^q)∂_qΨ)
"""
function Tj!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, B::BFieldT, met::MetT, Γ::Array{Float64, 2}, dΓ::Array{Float64, 3}, K::Array{Float64})

    #scale factor for this term.
    sf = - met.J[1] * jparonB(B, met) / 2.0

    #computes the Γ matrix and its derivative dΓ
    compute_Γ!(B, met, Γ, dΓ)


    for n in 1:3

        K .= 0.0

        #computes K, contribution to W for ∂^2 terms without Γ derivatives.
        compute_K!(met, Γ, K, n)

        #this loop computes the contributions for first derivatives,
        #so contains the derivatives of Γ.
        for q in 1:3
            #stores the intermediate example for each q.
            val = 0
            for i in 1:3, j in 1:3, k in 1:3
                val += Γ[i, n] / met.J[1] * lct[i, j, k] * dΓ[k, q, j]
            end

            #transpose is added to reflect that the two terms of Tj are identical except Ψ -> Φ.
            W[n, q] += val * sf
            W[q, n] += val * sf
        end

        #add the K contributions.
        #transpose is added to reflect that the two terms of Tj are identical except Ψ -> Φ.
        W[n, 4:9] .+= K .* sf
        W[4:9, n] .+= K .* sf

    end

end


"""
    function compute_K!(met::MetT, Γ::Array{Float64}, K::Array{Float64}, n::Int64)

Computes the K vector, which stores the contribution for Tj for the double derivative terms, i.e. W[4:9]
"""
function compute_K!(met::MetT, Γ::Array{Float64}, K::Array{Float64}, n::Int64)

    #double terms occur due to the symmetry in double derivatives
    #single terms when derivatives are with respect to the same variable twice.
    for i in 1:3, k in 1:3
        #drdr
        K[1] += Γ[i, n] * lct[i, 1, k] * Γ[k, 1] / met.J[1] 

        #drdθ, dθdr
        K[2] += Γ[i, n] * (lct[i, 1, k] * Γ[k, 2] + lct[i, 2, k] * Γ[k, 1]) / met.J[1]

        #drdζ, dζdr
        K[3] += Γ[i, n] * (lct[i, 1, k] * Γ[k, 3] + lct[i, 3, k] * Γ[k, 1]) / met.J[1]

        #dθdθ
        K[4] += Γ[i, n] * lct[i, 2, k] * Γ[k, 2] / met.J[1]

        #dθdζ, dζdθ
        K[5] += Γ[i, n] * (lct[i, 2, k] * Γ[k, 3] + lct[i, 3, k] * Γ[k, 2]) / met.J[1]

        #dζdζ
        K[6] += Γ[i, n] * lct[i, 3, k] * Γ[k, 3] / met.J[1]
    end
end


"""
    function compute_Γ!(B::BFieldT, met::MetT, Γ::Array{Float64, 2}, dΓ::Array{Float64, 3})

Computes the Γ matrix and its derivative dΓ. These are short hands for expressions appearing in Tj.
"""
function compute_Γ!(B::BFieldT, met::MetT, Γ::Array{Float64, 2}, dΓ::Array{Float64, 3})
    #may be able to compute Γ with gl * D...

    #Γ_i^j = δ^j_i - g_ik b^k b^j
    #dΓ_ik^j = ∂_k (Γ_i^j) = -b^l b^j ∂_k(g_il) - g_il b^l ∂_k(b^j) - b^j g_il ∂_k(b^l)
    #l contraction has been replaced with dot.
    for j in 1:3, i in 1:3
        Γ[i, j] = dδ(i, j) - dot(met.gl[i, :], B.b[:])*B.b[j]
        for k in 1:3
            dΓ[i, j, k] = (- B.b[j] * dot(B.b[:], met.dgl[i, :, k]) 
                        - dot(met.gl[i, :], B.b[:]) * B.db[j, k]
                        - B.b[j] * dot(met.gl[i, :], B.db[:, k]))

        end
    end
end


"""
    jparonB(B::BFieldT, met::MetT)

Computes the parrallel current divided by the magnitude of B, needed for the current term of W.
"""
function jparonB(B::BFieldT, met::MetT)
  
    jpar = 0.0

    #J_∥ = g_{ij}B^i J^j
    #J^j = (∇×B)^j = 1/J * ϵ^{jkl}∂_k B_l = 1/J * ϵ^{jkl}∂_k (g_{lp} B^p)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3, m in 1:3
        jpar += met.gl[i, j] * B.b[i] * 1.0 / met.J[1] * lct[j, k, l] * (met.gl[l, m] * B.dB[m, k] + met.dgl[l, m, k] * B.B[m])
    end
    return jpar/B.mag_B[1]
end


"""
    function get_lc_tensor()
        
Function for the levi-civita tensor.
"""
function get_lc_tensor()
    #probably a more elegant way to do this.
    lc = zeros(Float64, 3, 3, 3)
    lc[1, 2, 3] = 1.0
    lc[2, 3, 1] = 1.0
    lc[3, 1, 2] = 1.0
    lc[3, 2, 1] = -1.0
    lc[1, 3, 2] = -1.0
    lc[2, 1, 3] = -1.0
    return lc
end


#define the levi-civita tensor for curls.
const lct = get_lc_tensor()



"""
    function dδ(i::Int64, j::Int64)

Function representing the dirac delta.
"""
function dδ(i::Int64, j::Int64)
    if i==j
        return 1.0
    else
        return 0.0
    end
end

