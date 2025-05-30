"""
Struct storing the quantities needed to transform the toroidal metric and magnetic field into the qfm form.

### Fields
- coords::Array{Float64, 1} - array of (s, ϑ, φ), stored as array to keep struct immutable.
- JM::Array{Float64, 2} - Jacobian matrix, given by ∂x^i/∂x^μ, i.e. derivative matrix of the toroidal variables with respect to the qfm variables.
- JM_inv::Array{Float64, 2} - Inverse Jacobian matrix.
- dJM::Array{Float64, 3} - Derivative with respect to qfm variables of Jacobian matrix.
- dJM_inv::Array{Float64, 3} - Derivative with respect to qfm variables of inverse Jacobian matrix.
"""
struct CoordTransformT
    coords :: Array{Float64, 1} 
    JM :: Array{Float64, 2} 
    JM_inv :: Array{Float64, 2} 
    dJM :: Array{Float64, 3}
    dJM_inv :: Array{Float64, 3}
    function CoordTransformT()
        new(zeros(3), zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3))
    end
end



"""
    coord_transform!(s::Float64, ϑ::Float64, φ::Float64, CT::CoordTransformT, surf_itp::SurfaceITPT)

Function that computes the transformation from (s, ϑ, φ) (qfm coords) into standard toroidal coordinates (r, θ, ζ). This includes determining the Jacobian Matrix and associated derivatives, stored in the CT struct. This allows us to full convert required data between coordinates.
"""
function coord_transform!(s::Float64, ϑ::Float64, φ::Float64, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT)

    #the sd (surface data) struct is not very clear! probably needs a new name.
    pqMpol = surf_itp.M
    pqNtor = surf_itp.N
    dim1 = pqMpol+1
    dim2 = 2*pqNtor + 1

    #perhaps this should be precomputed, given we do it over and over again.
    mlist = collect(range(0, pqMpol))
    nlist = [collect(0:pqNtor);collect(-pqNtor:-1)] 

    #computes α
    for j in 1:1:length(nlist), i in 1:1:length(mlist)
        sd.α[i, j] = mlist[i] * ϑ - nlist[j] * φ
    end

    sd.cosα .= cos.(sd.α)
    sd.sinα .= sin.(sd.α)

    #this could arguably be made more efficient, as this is being called multiple times for the same s.
    #don't think this is a bottleneck though!
    itp_mat!(surf_itp, sd, s)
    
    #starting point of the coords, (r, θ, ζ)
    CT.coords .= [0.0, ϑ, φ]

    #starting point of JM matrix.
    #the jacobian Matrix, J_μ^i = ∂x^i/∂x^μ
    #JM[i, μ] denotes i as original var, (r, θ, ζ), and μ as new vars (s, ϑ, φ)
    #CT.JM .= [drds drdϑ drdφ; dθds dθdϑ dθdφ; dζds dζdϑ dζdφ]
    CT.JM .= [0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

    #starting point of dJM, all set to zero.
    #CT.dJM[:, :, 1] = [d2rdsds d2rdsdϑ d2rdsdφ; d2θdsds d2θdsdϑ d2θdsdφ; d2ζdsds d2ζdsdϑ d2ζdsdφ]
    #CT.dJM[:, :, 2] = [d2rdsdϑ d2rdϑdϑ d2rdϑdφ; d2θdsdϑ d2θdϑdϑ d2θdϑdφ; d2ζdsdϑ d2ζdϑdϑ d2ζdϑdφ]
    #CT.dJM[:, :, 3] = [d2rdsdφ d2rdϑdφ d2rdφdφ; d2θdsdφ d2θdϑdφ d2θdφdφ; d2ζdsdφ d2ζdϑdφ d2ζdφdφ]
    CT.dJM .= 0.0

    #perhaps less clear now, but maybe better?
    for j in 1:dim2, i in 1:dim1
        CT.coords[1] += sd.rcos[i, j] * sd.cosα[i, j] #r
        CT.coords[2] += sd.θsin[i, j] * sd.sinα[i, j] #θ

        CT.JM[1, 1] += sd.drcosds[i, j] * sd.cosα[i, j] #drds
        CT.JM[1, 2] += -mlist[i] * sd.rcos[i, j] * sd.sinα[i, j] #drdϑ
        CT.JM[1, 3] += nlist[j] * sd.rcos[i, j] * sd.sinα[i, j] #drdφ

        CT.JM[2, 1] += sd.dθsinds[i, j] * sd.sinα[i, j] #dθds
        CT.JM[2, 2] += mlist[i] * sd.θsin[i, j] * sd.cosα[i, j] #dθdϑ
        CT.JM[2, 3] += -nlist[j] * sd.θsin[i, j] * sd.cosα[i, j] #dθdφ

        CT.dJM[1, 1, 1] += sd.d2rcosdsds[i, j] * sd.cosα[i, j] #d2rdsds
        CT.dJM[1, 1, 2] += -mlist[i] * sd.drcosds[i, j] * sd.sinα[i, j] #d2rdsdϑ
        CT.dJM[1, 1, 3] += nlist[j] * sd.drcosds[i, j] * sd.sinα[i, j] #d2rdsdφ

        CT.dJM[1, 2, 2] += -mlist[i]^2 * sd.rcos[i, j] * sd.cosα[i, j] #d2rdϑdϑ
        CT.dJM[1, 2, 3] += mlist[i] * nlist[j] * sd.rcos[i, j] * sd.cosα[i, j] #d2rdϑdφ
        CT.dJM[1, 3, 3] += -nlist[j]^2 * sd.rcos[i, j] * sd.cosα[i, j] #d2rdφdφ

        CT.dJM[2, 1, 1] += sd.d2θsindsds[i, j] * sd.sinα[i, j] #d2θdsds
        CT.dJM[2, 1, 2] += mlist[i] * sd.dθsinds[i, j] * sd.cosα[i, j] #d2θdsdϑ
        CT.dJM[2, 1, 3] += -nlist[j] * sd.dθsinds[i, j] * sd.cosα[i, j] #d2θdsdφ

        CT.dJM[2, 2, 2] += -mlist[i]^2 * sd.θsin[i, j] * sd.sinα[i, j] #d2θdϑdϑ
        CT.dJM[2, 2, 3] += mlist[i] * nlist[j] * sd.θsin[i, j] * sd.sinα[i, j] #d2θdϑdφ
        CT.dJM[2, 3, 3] += -nlist[j]^2 * sd.θsin[i, j] * sd.sinα[i, j] #d2θdφdφ

    end

    #inverse
    CT.JM_inv .= inv(CT.JM)

    #symmetric second derivatives for dJM.
    CT.dJM[1, 2, 1] = CT.dJM[1, 1, 2] #d2rdsdϑ
    CT.dJM[1, 3, 1] = CT.dJM[1, 1, 3] #d2rdsdφ
    CT.dJM[2, 2, 1] = CT.dJM[2, 1, 2] #d2θdsdϑ
    CT.dJM[2, 3, 1] = CT.dJM[2, 1, 3] #d2θdsdφ
    CT.dJM[1, 3, 2] = CT.dJM[1, 2, 3] #d2rdϑdφ
    CT.dJM[2, 3, 2] = CT.dJM[2, 2, 3] #d2θdϑdφ

    #computes the derivative of the inverse transformation.
    #noteing that (K^-1)' = -(K^-1) @ K' @ (K^-1)
    for i in 1:3
        CT.dJM_inv[:, :, i] = - CT.JM_inv * CT.dJM[:, :, i] * CT.JM_inv
    end

end



"""
    B_transform!(tor_B::BFieldT, qfm_B::BFieldT, qfm_met::MetT, CT::CoordTransformT)

Function that fill out the BFieldT struct for qfm coordinates given the original BFieldT in toroidal coordinaets.
"""
function met_transform!(tor_met::MetT, qfm_met::MetT, CT::CoordTransformT)

   
    #reset the matrices.
    qfm_met.gl .= 0.0
    qfm_met.gu .= 0.0
    qfm_met.dgl .= 0.0

    #μ etc is used for indexing the new coords, i etc for old coords.
    #order of these loops can probably be improved.
    for i in 1:3, j in 1:3, μ in 1:3, ν in 1:3
        #JM[i, μ] = J^i_μ
        #g_{μν} = JM^i_μ g_{ij} JM^j_ν
        qfm_met.gl[μ, ν] += CT.JM[i, μ] * tor_met.gl[i, j] * CT.JM[j, ν] 

        #note that for inverse JM_inv[μ, i] = (J^i_μ)^{-1} = (J^{-1})^μ_i
        #explaination above is not clear, but we have to swap the indicies here...
        #index flipping here is for the inverse.
        #g^{μν} = JM_i^μ g_{ij} JM_j^ν
        qfm_met.gu[μ, ν] += CT.JM_inv[μ, i] * tor_met.gu[i, j] * CT.JM_inv[ν, j]

        for σ in 1:3

            #∂_σ g_{μν} = ∂_σ(JM^i_μ) g_{ij} JM^j_ν + JM^i_μ g_{ij} ∂_σ(JM^j_ν) 
            #+ JM^i_μ (JM^k_σ ∂_k(g_{ij})) JM^j_ν 
            #where the extra JM comes from the need to swap dg_tor to be a derivative in terms of hte new coordinates, not the old.
            qfm_met.dgl[μ, ν, σ] += (CT.dJM[i, μ, σ] * tor_met.gl[i, j] * CT.JM[j, ν]
                                     + CT.JM[i, μ] * tor_met.gl[i, j] * CT.dJM[j, ν, σ]
                                     + CT.JM[i, μ] * dot(CT.JM[:, σ], tor_met.dgl[i, j, :]) * CT.JM[j, ν])
                            
        end
    end
    #using derivtive of inverse matrix
    #(K^-1)' = - K^-1 (K)' K^-1
    for i in 1:3
        qfm_met.dgu[:, :, i] =  -1 .* qfm_met.gu * qfm_met.dgl[:, :, i] * qfm_met.gu
    end
    #perhap a try catche here?
    #sqrt here is annoying for values v close to zero
    #display((CT.coords[1], CT.coords[2], CT.coords[3]))
    qfm_met.J[1] = sqrt(det(qfm_met.gl))
    #qfm_met.J[1] = sqrt(abs(det(qfm_met.gl)))

    # now we take the derivative, noting that 
    # ∂(det(A)) = det(A) * Tr(A^{-1} ∂(A))
    # Here we just iterate over the derivatives, no funny buisness, as this is all in terms of the new metric.
    for μ in 1:3
        qfm_met.dJ[μ] = qfm_met.J[1]/2 * tr(qfm_met.gu * qfm_met.dgl[:, :, μ])
    end

end


"""
    B_transform!(tor_B::BFieldT, qfm_B::BFieldT, qfm_met::MetT, CT::CoordTransformT)

Function that fill out the BFieldT struct for qfm coordinates given the original BFieldT in toroidal coordinaets.
"""
function B_transform!(tor_B::BFieldT, qfm_B::BFieldT, qfm_met::MetT, CT::CoordTransformT)

    #B^μ = JM^μ_i B^i
    #Note that JM^μ_i is the inverse
    qfm_B.B .= CT.JM_inv * tor_B.B

    qfm_B.dB .= 0.0

    for σ in 1:3, μ in 1:3, i in 1:3
        #∂_σ B^μ = (∂_σ JM^μ_i) B^i + JM^μ_i * (JM^k_σ ∂_k B^i)
        #Extra JM (not inverse!) comes from converting derivative to new coords
        #dot is across k.
        qfm_B.dB[μ, σ] += CT.dJM_inv[μ, i, σ] * tor_B.B[i] + CT.JM_inv[μ, i] * dot(CT.JM[:, σ], tor_B.dB[i, :])
    end

    #same as normal, all in terms of new B and new met.
    Equilibrium.magnitude_B!(qfm_B, qfm_met)
    for i in 1:3
        qfm_B.b[i] = qfm_B.B[i]/qfm_B.mag_B[1]
    end
end

