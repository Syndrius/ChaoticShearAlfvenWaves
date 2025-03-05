
#Testing coordinate transform stuff with simpler example
#of cylindrical and spherical.
using LinearAlgebra
using MID

#not impossible that this will need to be bigger,
#i.e. we may want to transform more coords at once??
#but that is kind of against our method for MID.
struct CoordTsfmT
    JM :: Array{Float64, 2} #∂x^μ/∂x^i, i.e. deriv matrix of new vars w.r.t old vars.
    JM_inv :: Array{Float64, 2}
    dJM :: Array{Float64, 3}
    function CoordTsfmT()
        new(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3))
    end
end


function B_transform(tor_B, chaos_B, CT)
    #going to be v hard to verify all the derivs and stuff.

    #so we just need to transform B and dB, the other values can be computed, as long as we have the metrici information!

    #v simple transformation!
    chaos_B.B = CT.JM * tor_B.B

    for i in 1:3

        chaos_B.dB[:, i] = CT.JM * tor_B.dB[:, i] + CT.dJM[:, :, i] * tor_B.B
    end

    magnitude_B!(chaos_B, chaos_met)
    #perhaps qfm_B is a better name.
    chaos_B.b[1] = chaos_B.B[1]/chaos_B.mag_B
    chaos_B.b[2] = chaos_B.B[2]/chaos_B.mag_B
    chaos_B.b[3] = chaos_B.B[3]/chaos_B.mag_B
end


#this function computes the coordinate transform from toroidal like coordinates to qfm surfaces.
#returns the new coord and a struct containing the transoformation invormation.
#this will need more infor about the qfm surgaces etc.
function coord_transform(r, θ, ζ)

    #Obvs won't be creating a new version of this everytime...
    CT = CoordTsfmT()
    #will have to generalise to lists etc.
    m = 3
    n = 2

    #this will be a function of surface interpolation.
    s1 = 0.4
    ϑ1 = 0.3

    #obvs made these up!
    ds1dr = 0.1
    dϑ1dr = 0.3
    d2s1drdr = 0.05
    d2ϑ1drdr = 0.01


    s = s1 * cos(m * θ - n * ζ)

    ϑ = ϑ1 * sin(m * θ - n * ζ)

    φ = ζ


    dsdr = ds1dr * cos(m * θ - n * ζ)
    display(dsdr)
    dsdθ = -s1 * m * sin(m * θ - n * ζ)
    dsdζ = s1 * n * sin(m * θ - n * ζ)

    dϑdr = dϑ1dr * sin(m * θ - n * ζ)
    dϑdθ = ϑ1 * m * cos(m * θ - n * ζ)
    dϑdζ = -ϑ1 * n * cos(m * θ - n * ζ)

    dφdr = 0.0
    dφdθ = 0.0
    dφdζ = 1.0

    d2sdrdr = d2s1drdr * cos(m * θ - n * ζ)
    d2sdrdθ = -ds1dr * m * sin(m * θ - n * ζ)
    d2sdrdζ = ds1dr * n * sin(m * θ - n * ζ)

    d2sdθdθ = -s1 * m^2 * cos(m * θ - n * ζ)
    d2sdθdζ = s1 * m * n * cos(m * θ - n * ζ)
    d2sdζdζ = -s1 * n^2 * cos(m * θ - n * ζ)


    d2ϑdrdr = d2ϑ1drdr * sin(m * θ - n * ζ)
    d2ϑdrdθ = dϑ1dr * m * cos(m * θ - n * ζ)
    d2ϑdrdζ = -dϑ1dr * n * cos(m * θ - n * ζ)

    d2ϑdθdθ = -ϑ1 * m^2 * sin(m * θ - n * ζ)
    d2ϑdθdζ = ϑ1 * m * n * sin(m * θ - n * ζ)
    d2ϑdζdζ = -ϑ1 * n^2 * sin(m * θ - n * ζ)

    d2φdrdr = 0.0
    d2φdrdθ = 0.0
    d2φdrdζ = 0.0

    d2φdθdθ = 0.0
    d2φdθdζ = 0.0
    d2φdζdζ = 0.0

    #may want a JM struct tbh. And perhasp a function that fills it in.
    #probably the nicest way to deal with this.
    CT.JM = [dsdr dsdθ dsdζ; dϑdr dϑdθ dϑdζ; dφdr dφdθ dφdζ]
    CT.JM_inv = inv(CT.JM)

    
    #perhaps not the most efficient, but this is clearest I think.
    CT.dJM[:, :, 1] = [d2sdrdr d2sdrdθ d2sdrdζ; d2ϑdrdr d2ϑdrdθ d2ϑdrdζ; d2φdrdr d2φdrdθ d2φdrdζ] 

    CT.dJM[:, :, 2] = [d2sdrdθ d2sdθdθ d2sdθdζ; d2ϑdrdθ d2ϑdθdθ d2ϑdθdζ; d2φdrdθ d2φdθdθ d2φdθdζ]

    CT.dJM[:, :, 3] = [d2sdrdζ d2sdθdζ d2sdζdζ; d2ϑdrdζ d2ϑdθdζ d2ϑdζdζ; d2φdrdζ d2φdθdζ d2φdζdζ]
 
    return s, ϑ, φ, CT

end

#perhaps this should be called geometry transform or something??
function tor_to_chaos(r, θ, ζ, CT)

    #simplified version to double check with analytical case
    #note we will be using the cylindrical limit of torus for simplicity.

    tor_met = MID.MetT()
    chaos_met = MID.MetT()

    MID.cylindrical_metric!(tor_met, r, θ, ζ, 1.0)

   
    #so pretty sure we do need this!
    
    for μ in 1:3, ν in 1:3, i in 1:3, j in 1:3
        #damn so indexing this will be a nightmare.
        #although I guess we just have to know which index is contracting with the metric.
        #then this should be clearer.
        #We probbally need to be consistent with the upper or lower first thing.
        #going by this, we have defined J to have upper first. Seems ok.
        chaos_met.gl[μ, ν] += CT.JM_inv[i, μ] * tor_met.gl[i, j] * CT.JM_inv[j, ν] 

        chaos_met.gu[μ, ν] += CT.JM[μ, i] * tor_met.gu[i, j] * CT.JM[ν, j]

        #because of this call we need to get JM_inv
        for σ in 1:3, k in 1:3

            chaos_met.dgu[μ, ν, σ] += CT.JM_inv[k, σ] * (
                            tor_met.dgu[i, j, k] * CT.JM[μ, i] * CT.JM[ν, j]
                            + CT.dJM[μ, i, k] * tor_met.gu[i, j] * CT.JM[ν, j]
                            + CT.JM[μ, i] * tor_met.gu[i, j] * CT.dJM[ν, j, k])
            
        end

    end

    #using derivtive of inverse matrix
    #(K^-1)' = - K^-1 (K)' K^-1

    chaos_met.dgl[:, :, 1] =  -1 .* chaos_met.gl * chaos_met.dgu[:, :, 1] * chaos_met.gl
    #perhap a loop?
    chaos_met.dgl[:, :, 2] = -1 .* chaos_met.gl * chaos_met.dgu[:, :, 2] * chaos_met.gl
    chaos_met.dgl[:, :, 3] = -1 .* chaos_met.gl * chaos_met.dgu[:, :, 3] * chaos_met.gl

    chaos_met.J = sqrt(det(chaos_met.gl))
    #this is a properly crazy formula, hard to tell if it is actually working properly as dJ_cyl = [1.0, 0.0, 0.0], which we get. but not very conclusive that it will work in a more complicated case.
    for μ in 1:3
        chaos_met.dJ[μ] = chaos_met.J/2 * tr(chaos_met.gu * chaos_met.dgl[:, :, μ])
    end

    #display(JM)
    return chaos_met

end


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
