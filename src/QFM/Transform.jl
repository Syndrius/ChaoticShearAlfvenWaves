
#less of a disaster, but still danger zone.

#perhaps a coord inner struct would be better?
struct CoordTsfmT
    coords :: Array{Float64, 1} #array of (s, ϑ, φ) #stored like this so we can change them!
    JM :: Array{Float64, 2} #∂x^μ/∂x^i, i.e. deriv matrix of new vars w.r.t old vars.
    JM_inv :: Array{Float64, 2}
    dJM :: Array{Float64, 3}
    function CoordTsfmT()
        new(zeros(3), zeros(3, 3), zeros(3, 3), zeros(3, 3, 3))
    end
end


function B_transform!(tor_B, qfm_B, qfm_met, CT)
    #going to be v hard to verify all the derivs and stuff.

    #so we just need to transform B and dB, the other values can be computed, as long as we have the metrici information!

    #v simple transformation!
    qfm_B.B = CT.JM * tor_B.B

    for i in 1:3

        qfm_B.dB[:, i] = CT.JM * tor_B.dB[:, i] + CT.dJM[:, :, i] * tor_B.B
    end

    MagneticField.magnitude_B!(qfm_B, qfm_met)
    for i in 1:3
        qfm_B.b[i] = qfm_B.B[i]/qfm_B.mag_B
    end
    #perhaps qfm_B is a better name.
    #chaos_B.b[1] = chaos_B.B[1]/chaos_B.mag_B
    #chaos_B.b[2] = chaos_B.B[2]/chaos_B.mag_B
    #chaos_B.b[3] = chaos_B.B[3]/chaos_B.mag_B
end


#this function computes the coordinate transform from toroidal like coordinates to qfm surfaces.
#returns the new coord and a struct containing the transoformation invormation.
#this will need more infor about the qfm surgaces etc.
function coord_transform!(r, θ, ζ, CT, surf)


    #will have to generalise to lists etc.
    m = 3
    n = 2

    pqMpol = 24
    pqNtor = 8
    dim1 = pqMpol+1
    dim2 = 2*pqNtor + 1
    mlist = collect(range(0, pqMpol))

    collect(-pqNtor:0)
    nlist = [collect(0:pqNtor);collect(-pqNtor:-1)] 

    α = zeros((length(mlist), length(nlist)))

    for i in 1:1:length(mlist)
        for j in 1:1:length(nlist)
            α[i, j] = mlist[i] * θ - nlist[j] * ζ
        end
    end

    cosα = cos.(α)
    sinα = sin.(α)

    #note scos and ϑsin here is becuase we may later
    #need to use ssin and ϑcos as well.
    scos, ϑsin = itp_mat(surf, r, deriv=0)
    #unsure if these shuold be ds or dr.
    dscosdr, dϑsindr = itp_mat(surf, r, deriv=1)
    d2scosdrdr, d2ϑsindrdr = itp_mat(surf, r, deriv=2)
    

    #this will be a function of surface interpolation.
    #s1 = 0.4
    #ϑ1 = 0.3

    #obvs made these up!
    #ds1dr = 0.1
    #dϑ1dr = 0.3
    #d2s1drdr = 0.05
    #d2ϑ1drdr = 0.01
    s = 0.0
    ϑ = 0.0
    φ = ζ #no transformation for this.

    #all derivative terms we need.
    #perhaps there is a better way to do this.
    #may be clearer to store this as an array, 
    #i.e ds = zeros(3)
    # dds = zeros(6)?
    #unsure, both are bad I think.
    dsdr = 0.0
    dsdθ = 0.0
    dsdζ = 0.0

    dϑdr = 0.0
    dϑdθ = 0.0
    dϑdζ = 0.0

    dφdr = 0.0
    dφdθ = 0.0
    dφdζ = 1.0


    d2sdrdr = 0.0
    d2sdrdθ = 0.0
    d2sdrdζ = 0.0

    d2sdθdθ = 0.0
    d2sdθdζ = 0.0
    d2sdζdζ = 0.0


    d2ϑdrdr = 0.0
    d2ϑdrdθ = 0.0
    d2ϑdrdζ = 0.0

    d2ϑdθdθ = 0.0
    d2ϑdθdζ = 0.0
    d2ϑdζdζ = 0.0

    d2φdrdr = 0.0
    d2φdrdθ = 0.0
    d2φdrdζ = 0.0

    d2φdθdθ = 0.0
    d2φdθdζ = 0.0
    d2φdζdζ = 0.0

    for i in 1:dim1
        for j in dim2
            s += scos[i, j] * cosα[i, j]
            ϑ += ϑsin[i, j] * sinα[i, j]

            dsdr += dscosdr[i, j] * cosα[i, j]
            dsdθ += -mlist[i] * scos[i, j] * sinα[i, j]
            dsdζ += nlist[j] * scos[i, j] * sinα[i, j]

            dϑdr += dϑsindr[i, j] * sinα[i, j]
            dϑdθ += mlist[i] * ϑsin[i, j] * cosα[i, j]
            dϑdζ += -nlist[j] * ϑsin[i, j] * cosα[i, j]


            d2sdrdr += d2scosdrdr[i, j] * cosα[i, j]
            d2sdrdθ += -mlist[i] * dscosdr[i, j] * sinα[i, j]
            d2sdrdζ += nlist[j] * dscosdr[i, j] * sinα[i, j]

            d2sdθdθ += -mlist[i]^2 * scos[i, j] * cosα[i, j]
            d2sdθdζ += mlist[i] * nlist[j] * scos[i, j] * cosα[i, j]
            d2sdζdζ += -nlist[j]^2 * scos[i, j] * cosα[i, j]


            d2ϑdrdr += dϑsindr[i, j] * sinα[i, j]
            d2ϑdrdθ += mlist[i] * dϑsindr[i, j] * cosα[i, j]
            d2ϑdrdζ += -nlist[j] * dϑsindr[i, j] * cosα[i, j]

            d2ϑdθdθ += -mlist[i]^2 * ϑsin[i, j] * sinα[i, j]
            d2ϑdθdζ += mlist[i] * nlist[j] * ϑsin[i, j] * sinα[i, j]
            d2ϑdζdζ += -nlist[j]^2 * ϑsin[i, j] * sinα[i, j]

        end
    end

    #display((r, θ, ζ))
    #may want a JM struct tbh. And perhasp a function that fills it in.
    #probably the nicest way to deal with this.
    CT.JM[:, :] .= [dsdr dsdθ dsdζ; dϑdr dϑdθ dϑdζ; dφdr dφdθ dφdζ]
    #display(CT.JM)
    #so all the interpoaltions outside the domain are set to zero
    #causing this to give an error.
    #extrapolation etc is going to cause problemos.
    #display(CT.JM)
    CT.JM_inv .= inv(CT.JM)

    
    #perhaps not the most efficient, but this is clearest I think.
    CT.dJM[:, :, 1] = [d2sdrdr d2sdrdθ d2sdrdζ; d2ϑdrdr d2ϑdrdθ d2ϑdrdζ; d2φdrdr d2φdrdθ d2φdrdζ] 

    CT.dJM[:, :, 2] = [d2sdrdθ d2sdθdθ d2sdθdζ; d2ϑdrdθ d2ϑdθdθ d2ϑdθdζ; d2φdrdθ d2φdθdθ d2φdθdζ]

    CT.dJM[:, :, 3] = [d2sdrdζ d2sdθdζ d2sdζdζ; d2ϑdrdζ d2ϑdθdζ d2ϑdζdζ; d2φdrdζ d2φdθdζ d2φdζdζ]
 
    #return s, ϑ, φ, CT

end

#perhaps this should be called geometry transform or something??
function met_transform!(tor_met, qfm_met, CT)

    #simplified version to double check with analytical case
    #note we will be using the cylindrical limit of torus for simplicity.
   
    #so pretty sure we do need this!
    
    for μ in 1:3, ν in 1:3, i in 1:3, j in 1:3
        #damn so indexing this will be a nightmare.
        #although I guess we just have to know which index is contracting with the metric.
        #then this should be clearer.
        #We probbally need to be consistent with the upper or lower first thing.
        #going by this, we have defined J to have upper first. Seems ok.
        qfm_met.gl[μ, ν] += CT.JM_inv[i, μ] * tor_met.gl[i, j] * CT.JM_inv[j, ν] 

        qfm_met.gu[μ, ν] += CT.JM[μ, i] * tor_met.gu[i, j] * CT.JM[ν, j]

        #because of this call we need to get JM_inv
        for σ in 1:3, k in 1:3

            qfm_met.dgu[μ, ν, σ] += CT.JM_inv[k, σ] * (
                            tor_met.dgu[i, j, k] * CT.JM[μ, i] * CT.JM[ν, j]
                            + CT.dJM[μ, i, k] * tor_met.gu[i, j] * CT.JM[ν, j]
                            + CT.JM[μ, i] * tor_met.gu[i, j] * CT.dJM[ν, j, k])
            
        end

    end

    #using derivtive of inverse matrix
    #(K^-1)' = - K^-1 (K)' K^-1

    for i in 1:3
        qfm_met.dgl[:, :, i] =  -1 .* qfm_met.gl * qfm_met.dgu[:, :, i] * qfm_met.gl
    end
    #perhap a loop?
    #chaos_met.dgl[:, :, 2] = -1 .* chaos_met.gl * chaos_met.dgu[:, :, 2] * chaos_met.gl
    #chaos_met.dgl[:, :, 3] = -1 .* chaos_met.gl * chaos_met.dgu[:, :, 3] * chaos_met.gl

    #perhap a try catche here?
    qfm_met.J = sqrt(det(qfm_met.gl))
    #this is a properly crazy formula, hard to tell if it is actually working properly as dJ_cyl = [1.0, 0.0, 0.0], which we get. but not very conclusive that it will work in a more complicated case.
    for μ in 1:3
        qfm_met.dJ[μ] = qfm_met.J/2 * tr(qfm_met.gu * qfm_met.dgl[:, :, μ])
    end

    #display(JM)

end