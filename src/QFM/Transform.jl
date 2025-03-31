
#less of a disaster, but still danger zone.

#perhaps a coord inner struct would be better?
struct CoordTsfmT
    coords :: Array{Float64, 1} #array of (s, ϑ, φ) #stored like this so we can change them!
    JM :: Array{Float64, 2} #∂x^μ/∂x^i, i.e. deriv matrix of new vars w.r.t old vars.
    JM_inv :: Array{Float64, 2} 
    dJM :: Array{Float64, 3}
    jac :: Array{Float64, 1}
    function CoordTsfmT()
        new(zeros(3), zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(1))
    end
end


function B_transform!(tor_B, qfm_B, qfm_met, CT)
    #going to be v hard to verify all the derivs and stuff.

    #so we just need to transform B and dB, the other values can be computed, as long as we have the metrici information!

    #v simple transformation!
    #B^μ = JM^μ_i B^i
    #Note that JM^μ_i is the inverse
    qfm_B.B .= CT.JM_inv * tor_B.B

    #temporary solution!
    dJMinv = zeros(3, 3, 3)
    #derivative of inverse formula again!
    for i in 1:3
        dJMinv[:, :, i] = - CT.JM_inv * CT.dJM[:, :, i] * CT.JM_inv
    end

    qfm_B.dB .= 0.0

    for σ in 1:3, μ in 1:3, i in 1:3
        #∂_σ B^μ = (∂_σ JM^μ_i) B^i + JM^μ_i * (JM^k_σ ∂_k B^i)
        #Extra JM (not inverse!) comes from converting derivative to new coords
        #dot is across k.
        qfm_B.dB[μ, σ] += dJMinv[μ, i, σ] * tor_B.B[i] + CT.JM_inv[μ, i] * dot(CT.JM[:, σ], tor_B.dB[i, :])
    end

    #same as normal, all in terms of new B and new met!
    Equilibrium.magnitude_B!(qfm_B, qfm_met)
    for i in 1:3
        qfm_B.b[i] = qfm_B.B[i]/qfm_B.mag_B[1]
    end
end


function test_s(s, deriv)

    #test case to see if interpolation is problemo
    pqMpol = 24
    pqNtor = 8
    dim1 = pqMpol+1
    dim2 = 2*pqNtor + 1
    rcos = zeros(dim1, dim2)
    θsin = zeros(dim1, dim2)

    a = 0.001
    b = 0.0005

    if deriv == 0
        
        rcos[3, 2] = a*s^4
        θsin[3, 2] = b*(2*s + s^3)
    elseif deriv == 1
        rcos[3, 2] = a*4 * s^3
        θsin[3, 2] = b*(2 + 3*s^2)
    elseif deriv == 2
        rcos[3, 2] = a* 12 * s^2
        θsin[3, 2] = b * 6*s
    end
    #=
    if deriv == 0
        
        rcos .= a*s^4
        θsin .= b*(2*s + s^3)
    elseif deriv == 1
        rcos .= a*4 * s^3
        θsin .= b*(2 + 3*s^2)
    elseif deriv == 2
        rcos .= a* 12 * s^2
        θsin .= b * 6*s
    end
    =#
    return rcos, θsin
end

#this function computes the coordinate transform from toroidal like coordinates to qfm surfaces.
#returns the new coord and a struct containing the transoformation invormation.
#this will need more infor about the qfm surgaces etc.
##this needs to be changed, we have the coord transform flipped.
#function coord_transform!(r, θ, ζ, CT, surf)
#This function computes the coordinates transform from qfm coordinates, (s, ϑ, φ), to regular toroidal coordinates, (r, θ, ζ).
#It also returns the Jacobian matrix of the transformation so that the metric and magnetic field
#can also be transformed.
function coord_transform!(s, ϑ, φ, CT, surf)

    

    #will have to generalise to lists etc.
    #m = 3
    #n = 2
    #obvs needs to be passed in by the surface obj.
    pqMpol = 24
    pqNtor = 8
    dim1 = pqMpol+1
    dim2 = 2*pqNtor + 1
    mlist = collect(range(0, pqMpol))

    #collect(-pqNtor:0)
    nlist = [collect(0:pqNtor);collect(-pqNtor:-1)] 

    α = zeros((length(mlist), length(nlist)))

    for i in 1:1:length(mlist)
        for j in 1:1:length(nlist)
            α[i, j] = mlist[i] * ϑ - nlist[j] * φ
        end
    end

    cosα = cos.(α)
    sinα = sin.(α)

    #note scos and ϑsin here is becuase we may later
    #need to use ssin and ϑcos as well.
    #what is this actually giving us?
    #perhaps this is givnig us fourier components not 
    #actual θ, ζ values.
    rcos, θsin = itp_mat(surf, s, deriv=0)
    #unsure if these shuold be ds or dr.
    drcosds, dθsinds = itp_mat(surf, s, deriv=1)
    d2rcosdsds, d2θsindsds = itp_mat(surf, s, deriv=2)
    #=
    rcos, θsin = test_s(s, 0)
    drcosds, dθsinds = test_s(s, 1)
    d2rcosdsds, d2θsindsds = test_s(s, 2)
    =#
    

    #this will be a function of surface interpolation.
    #s1 = 0.4
    #ϑ1 = 0.3

    #obvs made these up!
    #ds1dr = 0.1
    #dϑ1dr = 0.3
    #d2s1drdr = 0.05
    #d2ϑ1drdr = 0.01
    r = 0.0
    #not 100% on this, but seems like it has to start here
    θ = ϑ #given the derivative term, this should probably start at ϑ or something?
    ζ = φ #no transformation for this.

    #all derivative terms we need.
    #perhaps there is a better way to do this.
    #may be clearer to store this as an array, 
    #i.e ds = zeros(3)
    # dds = zeros(6)?
    #unsure, both are bad I think.
    drds = 0.0
    drdϑ = 0.0 
    drdφ = 0.0

    dθds = 0.0
    dθdϑ = 1.0 #this starts at 1.0 because there must be an additional constant or something in the fourier expansion.
    dθdφ = 0.0

    dζds = 0.0
    dζdϑ = 0.0
    dζdφ = 1.0


    d2rdsds = 0.0
    d2rdsdϑ = 0.0
    d2rdsdφ = 0.0

    d2rdϑdϑ = 0.0
    d2rdϑdφ = 0.0
    d2rdφdφ = 0.0


    d2θdsds = 0.0
    d2θdsdϑ = 0.0
    d2θdsdφ = 0.0

    d2θdϑdϑ = 0.0
    d2θdϑdφ = 0.0
    d2θdφdφ = 0.0

    d2ζdsds = 0.0
    d2ζdsdϑ = 0.0
    d2ζdsdφ = 0.0

    d2ζdϑdϑ = 0.0
    d2ζdϑdφ = 0.0
    d2ζdφdφ = 0.0

    for i in 1:dim1
        for j in 1:dim2
            r += rcos[i, j] * cosα[i, j]
            θ += θsin[i, j] * sinα[i, j]

            drds += drcosds[i, j] * cosα[i, j]
            drdϑ += -mlist[i] * rcos[i, j] * sinα[i, j]
            drdφ += nlist[j] * rcos[i, j] * sinα[i, j]

            dθds += dθsinds[i, j] * sinα[i, j]
            dθdϑ += mlist[i] * θsin[i, j] * cosα[i, j]
            dθdφ += -nlist[j] * θsin[i, j] * cosα[i, j]


            d2rdsds += d2rcosdsds[i, j] * cosα[i, j]
            d2rdsdϑ += -mlist[i] * drcosds[i, j] * sinα[i, j]
            d2rdsdφ += nlist[j] * drcosds[i, j] * sinα[i, j]

            d2rdϑdϑ += -mlist[i]^2 * rcos[i, j] * cosα[i, j]
            d2rdϑdφ += mlist[i] * nlist[j] * rcos[i, j] * cosα[i, j]
            d2rdφdφ += -nlist[j]^2 * rcos[i, j] * cosα[i, j]


            d2θdsds += d2θsindsds[i, j] * sinα[i, j]
            d2θdsdϑ += mlist[i] * dθsinds[i, j] * cosα[i, j]
            d2θdsdφ += -nlist[j] * dθsinds[i, j] * cosα[i, j]

            d2θdϑdϑ += -mlist[i]^2 * θsin[i, j] * sinα[i, j]
            d2θdϑdφ += mlist[i] * nlist[j] * θsin[i, j] * sinα[i, j]
            d2θdφdφ += -nlist[j]^2 * θsin[i, j] * sinα[i, j]

        end
    end

    #th jacobian Matrix, J_μ^i = ∂x^i/∂x^μ
    #CT.JM .= [dsdr dsdθ dsdζ; dϑdr dϑdθ dϑdζ; dφdr dφdθ dφdζ]
    #JM[i, μ] denotes i as original var, (r, θ, ζ), and μ as new vars (s, ϑ, φ)
    CT.JM .= [drds drdϑ drdφ; dθds dθdϑ dθdφ; dζds dζdϑ dζdφ]
    #display(CT.JM)
    #so all the interpoaltions outside the domain are set to zero
    #causing this to give an error.
    #extrapolation etc is going to cause problemos.
    #display(CT.JM)
    ##don't actually think we need this anymore.
    #Don't think we ever actually needed this lol.
    CT.JM_inv .= inv(CT.JM)

    #this function is slow af (not surprising!) probably need to predefine arrays and reuse them!
    #or just add directly into these arrays? bit of a mess tbh.

    CT.dJM[:, :, 1] = [d2rdsds d2rdsdϑ d2rdsdφ; d2θdsds d2θdsdϑ d2θdsdφ; d2ζdsds d2ζdsdϑ d2ζdsdφ] 

    CT.dJM[:, :, 2] = [d2rdsdϑ d2rdϑdϑ d2rdϑdφ; d2θdsdϑ d2θdϑdϑ d2θdϑdφ; d2ζdsdϑ d2ζdϑdϑ d2ζdϑdφ]

    CT.dJM[:, :, 3] = [d2rdsdφ d2rdϑdφ d2rdφdφ; d2θdsdφ d2θdϑdφ d2θdφdφ; d2ζdsdφ d2ζdϑdφ d2ζdφdφ]

    #note that this makes assumptions on the form of the transformation,
    #namely that the ζ transform is simple
    CT.jac[1] = drds * dθdϑ - dθds * drdϑ
    #CT.jac[1] = sqrt(det(CT.JM))
    #perhaps we want include the derivative here?
    #
    CT.coords .= [r, θ, ζ]
    
    #perhaps not the most efficient, but this is clearest I think.
    #CT.dJM[:, :, 1] = [d2sdrdr d2sdrdθ d2sdrdζ; d2ϑdrdr d2ϑdrdθ d2ϑdrdζ; d2φdrdr d2φdrdθ d2φdrdζ] 

    #CT.dJM[:, :, 2] = [d2sdrdθ d2sdθdθ d2sdθdζ; d2ϑdrdθ d2ϑdθdθ d2ϑdθdζ; d2φdrdθ d2φdθdθ d2φdθdζ]

    #CT.dJM[:, :, 3] = [d2sdrdζ d2sdθdζ d2sdζdζ; d2ϑdrdζ d2ϑdθdζ d2ϑdζdζ; d2φdrdζ d2φdθdζ d2φdζdζ]
 
    #return s, ϑ, φ, CT

end

#perhaps this should be called geometry transform or something??
function met_transform!(tor_met, qfm_met, CT)

    #simplified version to double check with analytical case
    #note we will be using the cylindrical limit of torus for simplicity.
   
    #so pretty sure we do need this!

    #this is actually really, bad.
    #but required!
    qfm_met.gl .= 0.0
    qfm_met.gu .= 0.0
    qfm_met.dgl .= 0.0

    
    #μ etc is used for indexing the new coords, i etc for old coords.
    for μ in 1:3, ν in 1:3, i in 1:3, j in 1:3
        #damn so indexing this will be a nightmare.
        #although I guess we just have to know which index is contracting with the metric.
        #then this should be clearer.
        #We probbally need to be consistent with the upper or lower first thing.
        #going by this, we have defined J to have upper first. Seems ok.
        #JM[i, μ] = J^i_μ
        #g_{μν} = JM^i_μ g_{ij} JM^j_ν
        qfm_met.gl[μ, ν] += CT.JM[i, μ] * tor_met.gl[i, j] * CT.JM[j, ν] 

        #note that for inverse JM_inv[μ, i] = (J^i_μ)^{-1} = (J^{-1})^μ_i
        #explaination above is not clear, but we have to swap the indicies here...
        #index flipping here is for the inverse.
        #g^{μν} = JM_i^μ g_{ij} JM_j^ν
        qfm_met.gu[μ, ν] += CT.JM_inv[μ, i] * tor_met.gu[i, j] * CT.JM_inv[ν, j]

        #because of this call we need to get JM_inv
        #this seems to be a bit wrong unfort. Unsure tbh.
        for σ in 1:3

            #∂_σ g_{μν} = ∂_σ(JM^i_μ) g_{ij} JM^j_ν + JM^i_μ g_{ij} ∂_σ(JM^j_ν) 
            #+ JM^i_μ (JM^k_σ ∂_k(g_{ij})) JM^j_ν 
            #where the extra JM comes from the need to swap dg_tor to be a derivative in terms of hte new coordinates, not the old.
            qfm_met.dgl[μ, ν, σ] += (CT.dJM[i, μ, σ] * tor_met.gl[i, j] * CT.JM[j, ν]
                                     + CT.JM[i, μ] * tor_met.gl[i, j] * CT.dJM[j, ν, σ])

            #this term does not need be done separatly!
            qfm_met.dgl[μ, ν, σ] += CT.JM[i, μ] * dot(CT.JM[:, σ], tor_met.dgl[i, j, :]) * CT.JM[j, ν]
                            
        end

    end

    #using derivtive of inverse matrix
    #(K^-1)' = - K^-1 (K)' K^-1

    for i in 1:3
        qfm_met.dgu[:, :, i] =  -1 .* qfm_met.gu * qfm_met.dgl[:, :, i] * qfm_met.gu
    end
    #perhap a loop?
    #chaos_met.dgl[:, :, 2] = -1 .* chaos_met.gl * chaos_met.dgu[:, :, 2] * chaos_met.gl
    #chaos_met.dgl[:, :, 3] = -1 .* chaos_met.gl * chaos_met.dgu[:, :, 3] * chaos_met.gl

    #perhap a try catche here?
    #sqrt here is annoying for values v close to zero
    qfm_met.J[1] = sqrt(det(qfm_met.gl))
    #alternatively, we could use
    #looks like this expression is actually wrong. at least derivative computation becomes wrong
    #bit unfor as now we have to use the sqrt.
    #may be worth re-looking at this at some stage.
    #qfm_met.J[1] = tor_met.J[1] * CT.jac[1]
    #However, this requires storing CT.jac
    # now we take the derivative, noting that 
    # ∂(det(A)) = det(A) * Tr(A^{-1} ∂(A))
    # Here we just iterate over the derivatives, no funny buisness, as this is all in terms of the new metric.
    for μ in 1:3
        qfm_met.dJ[μ] = qfm_met.J[1]/2 * tr(qfm_met.gu * qfm_met.dgl[:, :, μ])
    end

    #display(JM)

end
