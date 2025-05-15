
#Testing coordinate transform stuff with simpler example
#of cylindrical and spherical.
using LinearAlgebra
using MID

#not impossible that this will need to be bigger,
#i.e. we may want to transform more coords at once??
#but that is kind of against our method for MID.
#%% Functions
struct CoordTsfmT
    coords :: Array{Float64, 1}
    JM :: Array{Float64, 2} #∂x^μ/∂x^i, i.e. deriv matrix of new vars w.r.t old vars.
    JM_inv :: Array{Float64, 2}
    dJM :: Array{Float64, 3}
    jac :: Array{Float64, 1} # unsure if this is actually useful at all.
    function CoordTsfmT()
        new(zeros(3), zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(1))
    end
end


#%%
#this function computes the coordinate transform from toroidal like coordinates to qfm surfaces.
#returns the new coord and a struct containing the transoformation invormation.
#this will need more infor about the qfm surgaces etc.
function coord_transform(s, ϑ, φ)

    #Obvs won't be creating a new version of this everytime...
    CT = CoordTsfmT()
    #will have to generalise to lists etc.
    m = 3
    n = 2

    #this will be a function of surface interpolation.
    r1 = 0.4
    θ1 = 0.3

    #obvs made these up!
    dr1ds = 0.1
    dθ1ds = 0.3
    d2r1dsds = 0.05
    d2θ1dsds = 0.01

    r = r1 * cos(m * ϑ - n * φ)

    θ = θ1 * sin(m * ϑ - n * φ)

    ζ = φ

    #first derivatives
    drds = dr1ds * cos(m * ϑ - n * φ)
    drdϑ = -r1 * m * sin(m * ϑ - n * φ)
    drdφ = r1 * n * sin(m * ϑ - n * φ)

    dθds = dθ1ds * sin(m * ϑ - n * φ)
    dθdϑ = θ1 * m * cos(m * ϑ - n * φ)
    dθdφ = -θ1 * n * cos(m * ϑ - n * φ)

    dζds = 0.0
    dζdϑ = 0.0
    dζdφ = 1.0

    #second derivatives of r
    d2rdsds = d2r1dsds * cos(m * ϑ - n * φ)
    d2rdsdϑ = -dr1ds * m * sin(m * ϑ - n * φ)
    d2rdsdφ = dr1ds * n * sin(m * ϑ - n * φ)

    d2rdϑdϑ = -r1 * m^2 * cos(m * ϑ - n * φ)
    d2rdϑdφ = r1 * m * n * cos(m * ϑ - n * φ)
    d2rdφdφ = -r1 * n^2 * cos(m * ϑ - n * φ)

    #second derivatives of θ
    d2θdsds = d2θ1dsds * sin(m * ϑ - n * φ)
    d2θdsdϑ = dθ1ds * m * cos(m * ϑ - n * φ)
    d2θdsdφ = -dθ1ds * n * cos(m * ϑ - n * φ)

    d2θdϑdϑ = -θ1 * m^2 * sin(m * ϑ - n * φ)
    d2θdϑdφ = θ1 * m * n * sin(m * ϑ - n * φ)
    d2θdφdφ = -θ1 * n^2 * sin(m * ϑ - n * φ)

    #second derivatives of ζ
    d2ζdsds = 0.0
    d2ζdsdϑ = 0.0
    d2ζdsdφ = 0.0

    d2ζdϑdϑ = 0.0
    d2ζdϑdφ = 0.0
    d2ζdφdφ = 0.0

    #may want a JM struct tbh. And perhasp a function that fills it in.
    #probably the nicest way to deal with this.
    CT.JM .= [drds drdϑ drdφ; dθds dθdϑ dθdφ; dζds dζdϑ dζdφ]
    CT.JM_inv .= inv(CT.JM)
    
    #doesn't really match Zhisong, but I think this is what he was going for.
    CT.jac[1] = drds * dθdϑ - dθds * drdϑ

    #perhaps not the most efficient, but this is clearest I think.
    CT.dJM[:, :, 1] = [d2rdsds d2rdsdϑ d2rdsdφ; d2θdsds d2θdsdϑ d2θdsdφ; d2ζdsds d2ζdsdϑ d2ζdsdφ] 

    CT.dJM[:, :, 2] = [d2rdsdϑ d2rdϑdϑ d2rdϑdφ; d2θdsdϑ d2θdϑdϑ d2θdϑdφ; d2ζdsdϑ d2ζdϑdϑ d2ζdϑdφ]

    CT.dJM[:, :, 3] = [d2rdsdφ d2rdϑdφ d2rdφdφ; d2θdsdφ d2θdϑdφ d2θdφdφ; d2ζdsdφ d2ζdϑdφ d2ζdφdφ]

    CT.coords .= [r, θ, ζ]
 
    return CT

end
#%%

#perhaps this should be called geometry transform or something??
function geometry_transform(s, ϑ, φ, CT)

    #simplified version to double check with analytical case
    #note we will be using the cylindrical limit of torus for simplicity.

    tor_met = MID.MetT()
    qfm_met = MID.MetT()

    #use the computed version of (r, θ, ζ) to compute the original metric.
    MID.cylindrical_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], 1.0)

    display(tor_met.gl)
    display(tor_met.dgl)

    #so pretty sure we do need this!
    #

    #note, JM[i, μ] = J^i_μ, makes index matching easier to think about!
    
    for μ in 1:3, ν in 1:3, i in 1:3, j in 1:3
        #damn so indexing this will be a nightmare.
        #although I guess we just have to know which index is contracting with the metric.
        #then this should be clearer.
        #We probbally need to be consistent with the upper or lower first thing.
        #going by this, we have defined J to have upper first. Seems ok.
        qfm_met.gl[μ, ν] += CT.JM[i, μ] * tor_met.gl[i, j] * CT.JM[j, ν] 

        #note that for inverse JM_inv[μ, i] = (J^i_μ)^{-1} = (J^{-1})^μ_i
        #explaination above is not clear, but we have to swap the indicies here...
        qfm_met.gu[μ, ν] += CT.JM_inv[μ, i] * tor_met.gu[i, j] * CT.JM_inv[ν, j]

        #because of this call we need to get JM_inv
        for σ in 1:3

            qfm_met.dgl[μ, ν, σ] += (CT.dJM[i, μ, σ] * tor_met.gl[i, j] * CT.JM[j, ν]
                                     + CT.JM[i, μ] * tor_met.gl[i, j] * CT.dJM[j, ν, σ])
#            for k in 1:3
            #could be combined now, this is the erivative of the met part, which has a bit more to it.
            qfm_met.dgl[μ, ν, σ] += dot(CT.JM[:, σ], tor_met.dgl[i, j, :]) * CT.JM[i, μ] * CT.JM[j, ν]
#            end

        end

    end

    for i in 1:3
        qfm_met.dgu[:, :, i] = -1 .* qfm_met.gu * qfm_met.dgl[:, :, i] * qfm_met.gu
    end

    #qfm_met.dgl[:, :, 1] = transpose(CT.dJM[:, :, 1])*tor_met.gl[:, :] * CT.JM + transpose(CT.JM) * tor_met.gl * CT.dJM[:, :, 1]
    
    #for μ in 1:3, ν in 1:3, σ in 1:3
        #order of CT.JM out the front is very important.
    #    qfm_met.dgl[μ, ν, σ] += dot(CT.JM[:, σ], tor_met.dgl[μ, ν, :])
    #end
    
    #for i in 1:3
    #    qfm_met.dgl[:, :, i] = transpose(CT.JM) * qfm_met.dgl[:, :, i] * CT.JM
    #end

    #using derivtive of inverse matrix
    #(K^-1)' = - K^-1 (K)' K^-1

    #qfm_met.dgu[:, :, 1] =  -1 .* qfm_met.gl * qfm_met.dgu[:, :, 1] * qfm_met.gl
    #perhap a loop?
    #qfm_met.dgu[:, :, 2] = -1 .* qfm_met.gl * qfm_met.dgu[:, :, 2] * qfm_met.gl
    #qfm_met.dgu[:, :, 3] = -1 .* qfm_met.gl * qfm_met.dgu[:, :, 3] * qfm_met.gl

    #think this is the easiest way to do this!
    #qfm_met.J = sqrt(det(qfm_met.gl))
    #unsure which is better, I guess the first way doesn't require we store CT.jac?
    #but this way we don't have to compute a sqrt(det())?
    qfm_met.J = tor_met.J * CT.jac[1]

    #this is a properly crazy formula, hard to tell if it is actually working properly as dJ_cyl = [1.0, 0.0, 0.0], which we get. but not very conclusive that it will work in a more complicated case.
    for μ in 1:3
        qfm_met.dJ[μ] = qfm_met.J/2 * tr(qfm_met.gu * qfm_met.dgl[:, :, μ])
        #qfm_met.dJ[μ] = tor_met.J * sqrt(det(CT.JM))
    end

    #display(JM)
    return qfm_met

end
#%%

#CT = CoordTsfmT()
CT = coord_transform(0.1, 0.8, 0.0)
display(CT.dJM)

#%%
qfm_met = geometry_transform(0.1, 0.8, 0.0, CT)
#display(qfm_met.gl)
#display(qfm_met.dgl)
#display(qfm_met.gu)
#display(qfm_met.dgu)
display(qfm_met.J)
display(qfm_met.dJ)

#%%


function B_transform(tor_B, chaos_B, CT, chaos_met)
    #going to be v hard to verify all the derivs and stuff.

    #so we just need to transform B and dB, the other values can be computed, as long as we have the metrici information!

    #v simple transformation!
    chaos_B.B = CT.JM_inv * tor_B.B

    dJMinv = zeros(3, 3, 3)

    #perhaps we start with this and see how we can avoid it later!
    dJMinv[:, :, 1] = -CT.JM_inv * CT.dJM[:, :, 1] * CT.JM_inv
    dJMinv[:, :, 2] = -CT.JM_inv * CT.dJM[:, :, 2] * CT.JM_inv
    dJMinv[:, :, 3] = -CT.JM_inv * CT.dJM[:, :, 3] * CT.JM_inv

    for σ in 1:3, μ in 1:3, i in 1:3
    #we may just have to back this one in, not really sure how best to test this.
    #better way to test maybe just to plot the poincare in the new coords, then we can see if it is straightish!
    #this is defs wrong though, need deriv w.r.t new var!
        chaos_B.dB[μ, σ] += dJMinv[μ, i, σ] * tor_B.B[i] + CT.JM_inv[μ, i] * dot(CT.JM[:, σ], tor_B.dB[i, :])
    end

    MID.MagneticField.magnitude_B!(chaos_B, chaos_met)
    #perhaps qfm_B is a better name.
    chaos_B.b[1] = chaos_B.B[1]/chaos_B.mag_B
    chaos_B.b[2] = chaos_B.B[2]/chaos_B.mag_B
    chaos_B.b[3] = chaos_B.B[3]/chaos_B.mag_B
end

tor_met = MID.MetT()
cylindrical_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], 1.0)
B_tor = MID.BFieldT()
#MID.compute_B!(B_tor, tor_met, chaos_q, MID.Structures.no_isl, MID.Structures.no_isl, CT.coords[1], CT.coords[2], CT.coords[3])
B_tor.B .= [CT.coords[1], CT.coords[2]/CT.coords[1], 1/ CT.coords[1]^2]
B_tor.dB[:, 1] = [1, -CT.coords[2]/CT.coords[1]^2, -2/ CT.coords[1]^3]
B_tor.dB[:, 2] = [0, 1/CT.coords[1], 0]


B_qfm = MID.BFieldT()
B_transform(B_tor, B_qfm, CT, qfm_met)
display(B_qfm.B)
display(B_qfm.dB[:, 1])
display(B_qfm.dB[:, 2])
display(B_qfm.dB[:, 3])

#%%


met = tor_to_chaos(0.1, 0.8, 0.0)

display(met.gu) #good!
display(met.gl) #good
display(met.J)
#not sure how to best verify the derivative terms
#however, they worked for simple sphr_to_cyl case, so hopefully ok
#also, for continuum we probably dont need them??
display(met.dJ) #doesn't match, hard to know if it is here or mathematica.
#same with this one, defs some similarities though!
#we can plausibly invert the transformation on mathematica
#then directly take the derivative of our new vars. unsure if mathematica will be able to do it thougH!
display(met.dgl[:, :, 1])

function sphr_to_cyl(r, θ, ζ)

    #(r, θ, ζ) -> (ρ, ϕ, z)

    ρ = r * sin(θ)
    ϕ = ζ
    z = r * cos(θ)

    dρdr = sin(θ)
    dρdθ = r*cos(θ)
    dρdζ = 0.0

    dϕdr = 0.0
    dϕdθ = 0.0
    dϕdζ = 1.0

    dzdr = cos(θ)
    dzdθ = -r*sin(θ)
    dzdζ = 0.0


    d2ρdrdr = 0.0
    d2ρdrdθ = cos(θ)
    d2ρdrdζ = 0.0

    d2ρdθdθ = -r*sin(θ)
    d2ρdθdζ = 0.0
    d2ρdζdζ = 0.0


    d2ϕdrdr = 0.0
    d2ϕdrdθ = 0.0
    d2ϕdrdζ = 0.0

    d2ϕdθdθ = 0.0
    d2ϕdθdζ = 0.0
    d2ϕdζdζ = 0.0


    d2zdrdr = 0.0
    d2zdrdθ = -sin(θ)
    d2zdrdζ = 0.0

    d2zdθdθ = -r*cos(θ)
    d2zdθdζ = 0.0
    d2zdζdζ = 0.0



    JM = [dρdr dρdθ dρdζ; dϕdr dϕdθ dϕdζ; dzdr dzdθ dzdζ]

    JM_inv = inv(JM)

    gl_sphr = [1.0 0.0 0.0 ; 0.0 r^2 0.0; 0.0 0.0 r^2*sin(θ)^2]

    gu_sphr = [1.0 0.0 0.0 ; 0.0 1/r^2 0.0; 0.0 0.0 1/(r^2*sin(θ)^2)]

    gl_cyl = zeros(3, 3)
    gu_cyl = zeros(3, 3)

    dgu_sphr = zeros(3, 3, 3)

    dgu_sphr[:, :, 1] = [0.0 0.0 0.0 ; 0.0 -2.0/r^3 0.0; 0.0 0.0 -2.0/(r^3*sin(θ)^2)]

    dgu_sphr[:, :, 2] = [0.0 0.0 0.0 ; 0.0 0.0 0.0; 0.0 0.0 -2.0 * cot(θ)/(r^2*sin(θ)^2)]

    dgu_cyl = zeros(3, 3, 3)

    dJM = zeros(3, 3, 3)

    #perhaps not the most efficient, but this is clearest I think.
    dJM[:, :, 1] = [d2ρdrdr d2ρdrdθ d2ρdrdζ; d2ϕdrdr d2ϕdrdθ d2ϕdrdζ; d2zdrdr d2zdrdθ d2zdrdζ] 

    dJM[:, :, 2] = [d2ρdrdθ d2ρdθdθ d2ρdθdζ; d2ϕdrdθ d2ϕdθdθ d2ϕdθdζ; d2zdrdθ d2zdθdθ d2zdθdζ]

    dJM[:, :, 3] = [d2ρdrdζ d2ρdθdζ d2ρdζdζ; d2ϕdrdζ d2ϕdθdζ d2ϕdζdζ; d2zdrdζ d2zdθdζ d2zdζdζ]

    J_cyl = 0.0
    J_sphr = r^2 * sin(θ)
    dJ_sphr = [2*r*sin(θ) r^2*cos(θ) 0.0]
    dJ_cyl = zeros(3)

    #don't think we actually want to do this!
    #need a proper expression for the det.
    #fortunaly one row is zeros!
    #this kind of ruins our general approach, but tbh we don't want to do all the derivs
    J_ct = -dϕdζ * (dρdr * dzdθ - dρdθ * dzdr)
    dJ_ct = zeros(3)

    #should be a product rule on first term but it is one so who cares.
    #how to we ensure the sign for this...
    dJ_ct[1] = -1.0 * (d2ρdrdr * dzdθ + dρdr * d2zdrdθ - (d2ρdrdθ * dzdr + dρdθ * d2zdrdr))

    dJ_ct[2] = -1.0 * (d2ρdrdθ * dzdθ + dρdr * d2zdθdθ - (d2ρdθdθ * dzdr + dρdθ * d2zdrdθ))
    dJ_ct[2] = -1.0 * (d2ρdrdζ * dzdθ + dρdr * d2zdθdζ - (d2ρdθdζ * dzdr + dρdθ * d2zdrdζ))
    #det(JM)

    #J_cyl = J_sphr * J_ct

    for μ in 1:3, i in 1:3
        dJ_cyl[μ] += JM_inv[i, μ] * (dJ_sphr[i] * J_ct + J_sphr * dJ_ct[i])
    end
    

    for μ in 1:3, ν in 1:3, i in 1:3, j in 1:3
        #damn so indexing this will be a nightmare.
        #although I guess we just have to know which index is contracting with the metric.
        #then this should be clearer.
        #We probbally need to be consistent with the upper or lower first thing.
        #going by this, we have defined J to have upper first. Seems ok.
        gl_cyl[μ, ν] += JM_inv[i, μ] * gl_sphr[i, j] * JM_inv[j, ν] 

        gu_cyl[μ, ν] += JM[μ, i] * gu_sphr[i, j] * JM[ν, j]

        #this is very rapidly getting out of hand.
        for σ in 1:3, k in 1:3

            dgu_cyl[μ, ν, σ] += JM_inv[k, σ] * (
                            dgu_sphr[i, j, k] * JM[μ, i] * JM[ν, j]
                            + dJM[μ, i, k] * gu_sphr[i, j] * JM[ν, j]
                            + JM[μ, i] * gu_sphr[i, j] * dJM[ν, j, k])
            #dgu_cyl[μ, ν, σ] +=  (
            #                dgu_sphr[i, j, k] * JM[μ, i] * JM[ν, j]
            #                + dJM[μ, i, k] * gu_sphr[i, j] * JM[ν, j]
            #                + JM[μ, i] * gu_sphr[i, j] * dJM[ν, j, k])
        end

    end

    #may need to take abs!
    
    gl_cyl2 = JM_inv' * gl_sphr * JM_inv

    #probs want to do this, then we don't have to invert JM.
    gl_cyl3 = inv(gu_cyl)
    gu_cyl2 = JM * gu_sphr * JM'

    dgl_cyl = zeros(3, 3, 3)

    #using derivtive of inverse matrix
    #(K^-1)' = - K^-1 (K)' K^-1

    dgl_cyl[:, :, 1] =  -1 .* gl_cyl * dgu_cyl[:, :, 1] * gl_cyl
    #perhap a loop?
    dgl_cyl[:, :, 2] = -1 .* gl_cyl * dgu_cyl[:, :, 2] * gl_cyl
    dgl_cyl[:, :, 3] = -1 .* gl_cyl * dgu_cyl[:, :, 3] * gl_cyl



    J_cyl = sqrt(det(gl_cyl))
    #this is a properly crazy formula, hard to tell if it is actually working properly as dJ_cyl = [1.0, 0.0, 0.0], which we get. but not very conclusive that it will work in a more complicated case.
    for μ in 1:3
        dJ_cyl[μ] = 1 /(2 * J_cyl) * J_cyl^2 * tr(gu_cyl * dgl_cyl[:, :, μ])
    end



    display("J_sphr is")
    display(J_sphr)

    display("dJ_sphr is")
    display(dJ_sphr)


    


    display("J_ct is")
    display(J_ct)

    display("dJ_ct is")
    display(dJ_ct)

    display("J_cyl is")
    display(J_cyl)

    display("dJ_cyl is")
    display(dJ_cyl)
    


    return ρ, ϕ, z

end


sphr_to_cyl(0.1, 0.8, 0.0)
