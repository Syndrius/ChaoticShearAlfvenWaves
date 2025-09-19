
"""
Struct for storing the metric which describes the geometry and related derivatives at each coordinate.

### Fields
- gl::Array{Float64, 2} - 3x3 matrix storing the lowered metric g_{ij}
- gu::Array{Float64, 2} - 3x3 matrix storing the raised metric g^{ij}
- dgl::Array{Float64, 3} - Derivative of gl, 3rd index labels coordinate that derivative is taken with respect to.
- dgl::Array{Float64, 3} - Derivative of gu, 3rd index labels coordinate that derivative is taken with respect to.
- J::Array{Float64} - Jacobian of the metric. Stored as an array so struct is immutable.
- dJ::Array{Float64, 1} - Derivative of J, index labels coordinate that derivative is taken with respect to.
"""
struct MetT
    gl :: Array{Float64, 2}
    gu :: Array{Float64, 2} 
    dgl :: Array{Float64, 3} 
    dgu :: Array{Float64, 3} 
    J :: Array{Float64, 1}
    dJ :: Array{Float64, 1} 
    function MetT()
        new(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3), zeros(1), zeros(3))
    end
end

#template metric function
function metric!(met::MetT, x1::Float64, x2::Float64, x3::Float64, R0::Float64)
end

"""
    toroidal_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)

Function that fills out the MetT struct for toroidal geometry. Metric elements taken from Energetic Particles in Tokamak Plasmas by Sergai Sharapov. Straight field line coordinates are radius (r), generalised poloidal angle (θ) and generalised toroidal angle (ζ), equal to negative of true toroidal angle. Additionally we assume low shear and approximate Δ' ≈ r/(4*R0).
"""
function rad_toroidal_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)
    
    #now lets try without the manual entries fk me. 
    Δp = r/(4*R0)
    Δpp = 1/(4*R0)

    ϵ = r/R0

    η = 1/2 * (ϵ+Δp)

    ηp = 1/2* (1/R0 + Δpp)

    met.J[1] = r * R0 * (1+2*ϵ*cos(θ))

    met.gl[1, 1] = 1-2*Δp * cos(θ)
    met.gl[1, 2] = r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gl[2, 1] = r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gl[2, 2] = r^2*(1+4*η*cos(θ) + 4*η^2)
    met.gl[3, 3] = R0^2*(1+2*ϵ*cos(θ))


    met.gu[:, :] = inv(met.gl)

    
    met.dJ[1] = R0 + 4*r * cos(θ)
    met.dJ[2] = -2 * r * R0*ϵ * sin(θ)

    #first two indicies give metric element, while third is derivative,
    #eg [1, 2, 3] is ∂g_{12}/∂ζ
    met.dgl[1, 1, 1] = -2*Δpp * cos(θ)
    met.dgl[1, 1, 2] = 2*Δp * sin(θ)

    met.dgl[1, 2, 1] = ((ϵ + Δp + r*Δpp) + r*(1/R0 + 2*Δpp)) * sin(θ)
    met.dgl[1, 2, 2] = r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgl[2, 1, 1] = ((ϵ + Δp + r*Δpp) + r*(1/R0 + 2*Δpp)) * sin(θ)
    met.dgl[2, 1, 2] = r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgl[2, 2, 1] = 2*r*(1+4*η*cos(θ) + 4*η^2) + r^2 * (4*ηp*cos(θ) + 8*η * ηp)
    met.dgl[2, 2, 2] = -r^2*(4*η*sin(θ))

    met.dgl[3, 3, 1] = 2*R0*cos(θ)
    met.dgl[3, 3, 2] = -2*R0^2*ϵ*sin(θ)

    for i in 1:3
        met.dgu[:, :, i] = - met.gu * met.dgl[:, :, i] * met.gu
    end

end


"""
    toroidal_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)

Function that fills out the MetT struct for toroidal geometry. Metric elements taken from Energetic Particles in Tokamak Plasmas by Sergai Sharapov. Straight field line coordinates are radius (r), generalised poloidal angle (θ) and generalised toroidal angle (ζ), equal to negative of true toroidal angle. Additionally we assume low shear and approximate Δ' ≈ r/(4*R0).
"""
function anal_toroidal_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)
    
    
    Δp = r/(4*R0)
    Δpp = 1/(4*R0)

    ϵ = r/R0

    η = 1/2 * (ϵ+Δp)

    ηp = 1/2* (1/R0 + Δpp)

    met.J[1] = r * R0 * (1+2*ϵ*cos(θ))

    met.gl[1, 1] = 1-2*Δp * cos(θ)
    met.gl[1, 2] = r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gl[2, 1] = r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gl[2, 2] = r^2*(1+4*η*cos(θ) + 4*η^2)
    met.gl[3, 3] = R0^2*(1+2*ϵ*cos(θ))


    met.gu[1, 1] = 1+2*Δp * cos(θ)
    met.gu[1, 2] = -1/r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gu[2, 1] = -1/r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gu[2, 2] = 1/r^2*(1-2*(ϵ+Δp)*cos(θ))
    met.gu[3, 3] = 1/R0^2*(1-2*ϵ*cos(θ))

    
    met.dJ[1] = R0 + 4*r * cos(θ)
    met.dJ[2] = -2 * r * R0*ϵ * sin(θ)

    #first two indicies give metric element, while third is derivative,
    #eg [1, 2, 3] is ∂g_{12}/∂ζ
    met.dgl[1, 1, 1] = -2*Δpp * cos(θ)
    met.dgl[1, 1, 2] = 2*Δp * sin(θ)

    met.dgl[1, 2, 1] = ((ϵ + Δp + r*Δpp) + r*(1/R0 + 2*Δpp)) * sin(θ)
    met.dgl[1, 2, 2] = r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgl[2, 1, 1] = ((ϵ + Δp + r*Δpp) + r*(1/R0 + 2*Δpp)) * sin(θ)
    met.dgl[2, 1, 2] = r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgl[2, 2, 1] = 2*r*(1+4*η*cos(θ) + 4*η^2) + r^2 * (4*ηp*cos(θ) + 8*η * ηp)
    met.dgl[2, 2, 2] = -r^2*(4*η*sin(θ))

    met.dgl[3, 3, 1] = 2*R0*cos(θ)
    met.dgl[3, 3, 2] = -2*R0^2*ϵ*sin(θ)

    met.dgu[1, 1, 1] = 2*Δpp * cos(θ)
    met.dgu[1, 1, 2] = -2*Δp * sin(θ)

    met.dgu[1, 2, 1] = (1/r^2 * (ϵ + Δp + r*Δpp) - 1/r*(1/R0 + 2*Δpp)) * sin(θ)
    met.dgu[1, 2, 2] = -1/r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgu[2, 1, 1] = (1/r^2 * (ϵ + Δp + r*Δpp) - 1/r*(1/R0 + 2*Δpp)) * sin(θ)
    met.dgu[2, 1, 2] = -1/r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgu[2, 2, 1] = -2/r^3 * (1-2*(ϵ+Δp)*cos(θ)) + 1/r^2 * (-2*(1/R0+Δpp)*cos(θ))
    met.dgu[2, 2, 2] = 2/r^2*(ϵ+Δp)*sin(θ)

    met.dgu[3, 3, 1] = -2*cos(θ)/R0^3
    met.dgu[3, 3, 2] = 2*ϵ*sin(θ)/R0^2


end

#starting to come together, looks to work for low m/n and not many κ values, but doesn't work
#for all cases.
function island_metric!(met::MetT, κ::Float64, ᾱ::Float64, τ::Float64, R0::Float64, isl::CoordIslandT)

    #first we compute our intermediate metric (κ, β, ζ)
    #where √κsin(β) = sin(m0 α/2)

    #we have, for the first transformation
    #κ = 4/w^2 * (ψ-ψ0)^2 + sin^2(m0*α/2) -> κ(ψ, θ, ζ)
    #β = arctan(w/2 * sin(m0*α/2)/(ψ-ψ0)) -> β(ψ, θ, ζ)
    #τ = ζ -> τ(ψ, θ, ζ)

    #We then straighten these coordinates via
    #ᾱ = π/(2*K(κ)) * F(β, κ) -> ᾱ(κ, θ, ζ)

    #other relations
    #ψ-ψ0 = w*sqrt(κ)/2 * cos(β)
    #but we will use β as an intermediate for derivatives to reduce complexity.

    #once we have the non-deriv metric
    #we should be able to replictate the island continuum
    #then we can check all the derivs with finite diff.



    #in first transformed coordinates
    #need to change to ab.
    #maybe even swap from ab to something else.
    #∇κ2 = 32*ψ*κ / w^2 * cos(β)^2 + bhat * κ * sin(β)^2

    #bhat = (m^2/(2ψ) + n^2 / R0^2) * (1-κ*sin(β)^2)

    #∇κ∇β = (-16ψ/w^2 + bhat/2) * cos(β) * sin(β)

    #∇κ∇τ = n/R0^2 * sqrt(κ) * sin(β) * sqrt(1-κ*sin^2(β))

    #∇β2 = 8*ψ*sin(β)^2 / (w^2*κ) + bhat/(2κ) * cos(β)^2

    #∇β∇τ = n/(2*sqrt(κ)*R0^2) * cos(β) * sqrt(1-κ*sin(β)^2) / R0^2 #this term is probbaly negligible.

    #∇τ2 = 1/R0^2


    K, E = Elliptic.ellipke(κ)

    arg = 2*K / π * ᾱ

    sn, cn, dn = Elliptic.ellipj(arg, κ)
    β = Elliptic.Jacobi.am(arg, κ)

    Eβ = Elliptic.E(β, κ)
    Fβ = Elliptic.F(β, κ)

    Z = Eβ - E / K * Fβ

    #extract the island info
    ψ0 = isl.ψ0
    w = isl.w
    m0 = isl.m0
    n0 = isl.n0

    ψ = ψ0 + w * sqrt(κ) / 2  * cn

    chat = dn #unsure if this is clearer!

    bhat = (m0^2/(2ψ) + n0^2 / R0^2) * chat^2

    ∇κ∇β = (-16ψ/w^2 + bhat/2) * cn * sn

    ∇β2 = 8*ψ*cn^2 / (w^2*κ) + bhat/(2κ) * cn^2

    ∇β∇τ = n0/(2*sqrt(κ)*R0^2) * cn * dn / R0^2 #this term is probbaly negligible.

    ∇κ2 = 32 * κ / w^2 * cn^2 + bhat * κ * sn^2

    dᾱdκ = π / (4*κ * (1-κ) * K) * (Z - κ*sn*cn / dn)

    dᾱdβ = π / (2 * K * chat)

    ∇κ∇ᾱ = dᾱdκ * ∇κ2 + dᾱdβ * ∇κ∇β

    ∇κ∇τ = n0/R0^2 * sqrt(κ) * sn * dn

    ∇ᾱ2 = (dᾱdκ)^2 * ∇κ2 + dᾱdκ * dᾱdβ * ∇κ∇β + (dᾱdβ)^2 * ∇β2

    ∇ᾱ∇τ = dᾱdκ * ∇κ∇τ + dᾱdβ * ∇β∇τ

    ∇τ2 = 1/ R0^2

    met.gu[1, 1] = ∇κ2
    met.gu[1, 2] = ∇κ∇ᾱ
    met.gu[1, 3] = ∇κ∇τ

    met.gu[2, 1] = ∇κ∇ᾱ
    met.gu[2, 2] = ∇ᾱ2
    met.gu[2, 3] = ∇ᾱ∇τ

    met.gu[3, 1] = ∇κ∇τ
    met.gu[3, 2] = ∇ᾱ∇τ
    met.gu[3, 3] = ∇τ2

    met.gl .= inv(met.gu)
    #display(sqrt(det(met.gl)))
    #display(w/(m0*π) * K * R0)
    #display((m0*π) / (w * K) * R0)
    #met.J[1] = sqrt(det(met.gl))
    met.J[1] = w/(m0*π) * K * R0


end 
#=

    #and some common vars used
    sκ = sqrt(κ)

    #we can then compute our original vars
    ψ = ψ0 + w/2 * sκ * cos(β)
    α = 2/m0 * asin(sκ*sin(β))

    r = sqrt(2ψ)

    #and metric elements, assuming cylinder.
    ∇ψ2 = r^2
    ∇θ2 = 1/r^2
    ∇ζ2 = 1/R0^2

    dκdψ = 8/w^2*(ψ-ψ0)
    display(dκdψ)
    dκdθ = m0*sin(m0*α/2)*cos(m0*α/2)
    dκdζ = n0*sin(m0*α/2)*cos(m0*α/2)

    display(dκdθ)
    display(dκdζ)
    dβdψ = -2/(w*κ^2) * sin(m0*α/2)
    dβdθ = m0 / (w*κ^2) * (ψ-ψ0) * cos(m0*α/2)
    dβdζ = n0 / (w*κ^2) * (ψ-ψ0) * cos(m0*α/2)

    #this probably shows it is possible to do the toroidal version, probbaly not very interesting.
    ∇κ2 = (dκdψ)^2 * ∇ψ2 + (dκdθ)^2 * ∇θ2 + (dκdζ)^2 * ∇ζ2
    #∇β2 = (dβdψ)^2 * ∇ψ2 + (dβdθ)^2 * ∇θ2 + (dβdζ)^2 * ∇ζ2
    #∇κ∇β = 

    #Axel's version of this is wildy more complicated, as he writes this in terms of ᾱ instead of β, may be simpler once we take derivatives?
    dᾱdκ = π / (8*(κ-1) * κ * K^2) * (2 * E * Fβ + K * (-2*Eβ + κ*sin(2*β) / sqrt(1-κ*sin(β)^2))) #last term can probably be simplified a bit.
    dᾱdβ = π / (2 * K * sqrt(1-κ*sin(β)^2))

    #we will need the derives of each of these, and all the terms that makes them up!
    dᾱdψ = dᾱdκ * dκdψ + dᾱdβ * dβdψ
    dᾱdθ = dᾱdκ * dκdθ + dᾱdβ * dβdθ
    dᾱdζ = dᾱdκ * dκdζ + dᾱdβ * dβdζ

    dτdζ = 1.0 #just kept for symmetry.

    #computing derivatives of this will be awful
    #perhaps we can follow a similar pattern with the inverse derivs
    #this should be the actual values of the metric tho! assuming no typos.
    ∇κ2 = (dκdψ)^2 * ∇ψ2 + (dκdθ)^2 * ∇θ2 + (dκdζ)^2 * ∇ζ2
    ∇κ∇ᾱ = dκdψ * dᾱdψ * ∇ψ2 + dκdθ * dᾱdθ * ∇θ2 + dκdζ * dᾱdζ * ∇ζ2
    ∇κ∇τ = dκdζ * dτdζ * ∇ζ2
    ∇ᾱ2 = (dᾱdψ)^2 * ∇ψ2 + (dᾱdθ)^2 * ∇θ2 + (dᾱdζ)^2 * ∇ζ2
    ∇ᾱ∇τ = dᾱdζ * dτdζ * ∇ζ2
    ∇τ2 = (dτdζ)^2 * ∇ζ2


    met.gu[1, 1] = ∇κ2
    met.gu[1, 2] = ∇κ∇ᾱ
    met.gu[1, 3] = ∇κ∇τ

    met.gu[2, 1] = ∇κ∇ᾱ
    met.gu[2, 2] = ∇ᾱ2
    met.gu[2, 3] = ∇ᾱ∇τ

    met.gu[3, 1] = ∇κ∇τ
    met.gu[3, 2] = ∇ᾱ∇τ
    met.gu[3, 3] = ∇τ2
end

    #this will be annoying though, as we will need to find ψ(κ, ᾱ, τ) etc.
    #unsure if this approach will be easier or harder than the one below
    #think this is better because we can check each step of the derivatives as we go with finite difference etc.
    d∇ψ2dκ = d∇ψ2dψ * dψdκ + d∇ψ2dθ * dθdκ

    #will this work?
    #plausible, can't see it working for ᾱ though!
    #may need to convert each of our expressions in terms of κ, ᾱ, τ...
    #nah think it should work, just be a mess
    d∇κ2dκ = 2*ddκdψdκ * ∇ψ2 + (dκdψ)^2 * d∇ψ2dκ







end
=#



"""
    island_metric!(met::MetT, κ::Float64, ᾱ::Float64, φ::Float64, R0::Float64, isl::IslandT)

Function that fills out the MetT structure for magnetic island striaght field line geometry. Based on Konies et al 2024. κ is the radial variable, measure from island center, ᾱ is the striaghtened helical angle and φ is a typical toroidal angle. This assumed the original geometry was cylindrical.
"""
function old_island_metric!(met::MetT, κ::Float64, ᾱ::Float64, φ::Float64, R0::Float64, isl::CoordIslandT)

    #derivs have all been checked
    #assuming original expressions are correct
    #then this is all good.
    #needs a big ol clean though obvs

    K, E = Elliptic.ellipke(κ)

    arg = 2*K / π * ᾱ

    sn, cn, dn = Elliptic.ellipj(arg, κ)


    β = Elliptic.Jacobi.am(arg, κ)

    Eβ = Elliptic.E(β, κ)
    Fβ = Elliptic.F(β, κ)

    Z = Eβ - E / K * Fβ



    #hopefully we don't actually need the island for this...
    #think we will, probably need an island submodule... big RIP for combining everything.
    #maybe we need an island coords flag or something to separate them all, or just an island grid, which is the same functionally,
    #but uses κ etc and has a metric function that accepts an island.
    r02 = (isl.r0)^2
    w = isl.w
    m0 = isl.m0
    n0 = isl.n0
    sqκ = sqrt(κ)
    dsqκ = 1/(2*sqκ)
    r2 = r02 + w * sqκ * cn




    #is this supposed to be something else?
    bhat = (m0^2/r2 + n0^2/R0^2) * dn^2

    #β is used the same as Axel, but we have replaced βs with ᾱ to match our eventual method.

    #this assumes cylindrical metric elements. Hence the r^2 and R0^2, for torus they are more complicated.
    ∇κ2 = 16*r2*κ*cn^2/w^2 + bhat*κ*sn^2

    #last bit comes from comparison with Axel
    ∇κ∇β = -8*r2/w^2 * cn*sn + bhat/2 * sn*cn 

    ∇β2 = 4*r2/(κ*w^2) * sn^2 + bhat/(4*κ) * cn^2

    #this is defs where the problemo is!!!
    #and or the definitions of elliptic arg.
    dᾱdκ =  π / (4*K * κ * (1-κ)) * (Z - κ*cn*sn / dn)
    
    dᾱdβ = π / (2*K * dn)

    cd = cn / dn
    sc = sn / cn


    dKdκ = (E - (1-κ)*K) / (2*(1-κ)*κ)
    dEdκ = (E - K)/(2κ)
    dargdκ = 2/π * ᾱ * dKdκ
    dargdᾱ = 2*K/π
    
    dsndκ = cn * dn * dargdκ - 1/(2*(1-κ)*κ)* cn * dn * (Eβ - (1-κ) * arg - κ * cd * sn)
    dsndᾱ = cn * dn * dargdᾱ

    dcndκ = 1/2 * dn * sn * (1/((1-κ)κ) * (Eβ - (1-κ)*arg - κ*cd*sn) -  2*dargdκ)
    dcndᾱ = -dn * sn * dargdᾱ

    ddndκ = 1/2 * cn * sn * (1/(1-κ) * (Eβ - (1-κ)*arg - dn*sc) -  2*κ*dargdκ)
    ddndᾱ = -κ*cn * sn * dargdᾱ

    dr2dκ = 1/2*w*κ^(-1/2)*cn + w*sqκ*dcndκ
    dr2dᾱ = w*sqκ*dcndᾱ

    dbhatdκ = -m0^2/r2^2 * dr2dκ * dn^2 + (m0^2/r2 + n0^2/R0^2) * 2 * dn * ddndκ
    dbhatdᾱ = -m0^2/r2^2 * dr2dᾱ * dn^2 + (m0^2/r2 + n0^2/R0^2) * 2 * dn * ddndᾱ

    
    met.gu[1, 1] = ∇κ2
    met.gu[1, 2] = dᾱdκ * ∇κ2 + dᾱdβ * ∇κ∇β
    met.gu[1, 3] = n0/R0^2 * sqκ * sn*dn

    met.gu[2, 1] = dᾱdκ * ∇κ2 + dᾱdβ * ∇κ∇β
    met.gu[2, 2] = dᾱdκ^2 * ∇κ2 + 2 * dᾱdβ * dᾱdκ * ∇κ∇β + dᾱdβ^2 * ∇β2
    met.gu[2, 3] = n0/R0^2*(dᾱdκ * sqκ * sn*dn + dᾱdβ * dsqκ * dn * cn)

    met.gu[3, 1] = n0/R0^2 * sqκ * sn*dn
    met.gu[3, 2] = n0/R0^2*(dᾱdκ * sqκ * sn*dn + dᾱdβ * dsqκ * dn * cn)
    met.gu[3, 3] = 1/R0^2

    d∇κ2dκ = (16/w^2 * (r2*cn^2 + dr2dκ*κ*cn^2 + r2*κ*2*cn*dcndκ) +
              (bhat*sn^2 + dbhatdκ*κ*sn^2 + bhat*κ*2*sn*dsndκ))
    d∇κ2dᾱ = (16/w^2 * κ * (dr2dᾱ*cn^2 + r2*2*cn*dcndᾱ) +
              κ * (dbhatdᾱ*sn^2 + bhat*2*sn*dsndᾱ))

    d∇β2dκ = (4/w^2*(dr2dκ/κ*sn^2 - r2/κ^2*sn^2 + 2*r2/κ*sn*dsndκ)+
              1/4*(dbhatdκ/κ*cn^2 - bhat/κ^2*cn^2 + 2*bhat/κ*cn*dcndκ))

    d∇β2dᾱ = (4/w^2*(dr2dᾱ/κ*sn^2 + 2*r2/κ*sn*dsndᾱ)+
              1/4*(dbhatdᾱ/κ*cn^2 + 2*bhat/κ*cn*dcndᾱ))

    dβdκ = - 1/(2*(1-κ)*κ) * ((Eβ - (1-κ)*arg)*dn - κ*cn * sn) + dn * dargdκ
    dβdᾱ = dn * dargdᾱ

    dEβdκ = (Eβ - Fβ) /(2κ) + sqrt(1-κ*sin(β)^2)*dβdκ
    dEβdᾱ = sqrt(1-κ*sin(β)^2)*dβdᾱ

    dFβdκ = 1/4*(2*Eβ/(κ-κ^2) - 2*Fβ/κ - sin(2β)/((1-κ)*sqrt(1-κ*sin(β)^2)) + 4*dβdκ/sqrt(1-κ*sin(β)^2))
    dFβdᾱ = dβdᾱ/sqrt(1-κ*sin(β)^2)

    dZdκ = dEβdκ - (dEdκ / K * Fβ - E * dKdκ * Fβ /K^2 + E/K * dFβdκ)
    dZdᾱ = dEβdᾱ - E/K * dFβdᾱ

    #fk me
    d2ᾱdκ2 = (π / (4*K * κ * (1-κ)) * (dZdκ - 1/dn*(cn*sn + κ*dcndκ*sn + κ*cn*dsndκ - κ*cn*sn*ddndκ/dn)) - 
              π * (Z - κ*cn*sn/dn) * (4*dKdκ*κ*(1-κ) + 4*K*(1-2κ)) / (4*K*κ*(1-κ))^2)

    #should this always be zero? this is an odd expression.
    d2ᾱdκdᾱ = π / (4*K * κ * (1-κ)) * (dZdᾱ - κ/dn*(dcndᾱ*sn + cn*dsndᾱ - cn*sn*ddndᾱ/dn))

    d2ᾱdβdκ = - π / (2*K*dn) * (dKdκ / K + ddndκ / dn)
    #perhaps this should just be zero as well?
    #this is not zero, guess because this is partial of a grad type expression, not partial of partial
    d2ᾱdβdᾱ = - π / (2*K*dn^2) * ddndᾱ

    d∇κ∇βdκ = (-8/w^2 * (dr2dκ*cn*sn + r2*dcndκ*sn + r2*cn*dsndκ) +
               1/2 * (dbhatdκ*cn*sn + bhat*dcndκ*sn + bhat*cn*dsndκ))

    d∇κ∇βdᾱ = (-8/w^2 * (dr2dᾱ*cn*sn + r2*dcndᾱ*sn + r2*cn*dsndᾱ) +
               1/2 * (dbhatdᾱ*cn*sn + bhat*dcndᾱ*sn + bhat*cn*dsndᾱ))


    met.dgu[1, 1, 1] = d∇κ2dκ
    met.dgu[1, 1, 2] = d∇κ2dᾱ

    met.dgu[1, 2, 1] = (d2ᾱdκ2 * ∇κ2 + dᾱdκ * d∇κ2dκ) + (d2ᾱdβdκ * ∇κ∇β + dᾱdβ * d∇κ∇βdκ)

    met.dgu[1, 2, 2] = (d2ᾱdκdᾱ * ∇κ2 + dᾱdκ * d∇κ2dᾱ) + (d2ᾱdβdᾱ * ∇κ∇β + dᾱdβ * d∇κ∇βdᾱ)

    #surprising that n0 is not squared tbh
    met.dgu[1, 3, 1] = n0/R0^2 * (dsqκ*sn*dn + sqκ*dsndκ*dn + sqκ*sn*ddndκ)
    met.dgu[1, 3, 2] = n0/R0^2 * (sqκ * dsndᾱ*dn + sqκ*sn*ddndᾱ)

    met.dgu[2, 1, 1] = (d2ᾱdκ2 * ∇κ2 + dᾱdκ * d∇κ2dκ) + (d2ᾱdβdκ * ∇κ∇β + dᾱdβ * d∇κ∇βdκ)

    met.dgu[2, 1, 2] = (d2ᾱdκdᾱ * ∇κ2 + dᾱdκ * d∇κ2dᾱ) + (d2ᾱdβdᾱ * ∇κ∇β + dᾱdβ * d∇κ∇βdᾱ)

    met.dgu[2, 2, 1] = (2*dᾱdκ*d2ᾱdκ2*∇κ2 + dᾱdκ^2*d∇κ2dκ + 
                        2*(d2ᾱdβdκ*dᾱdκ*∇κ∇β + dᾱdβ*d2ᾱdκ2*∇κ∇β + dᾱdβ*dᾱdκ*d∇κ∇βdκ) +
                        2*dᾱdβ*∇β2*d2ᾱdβdκ + dᾱdβ^2*d∇β2dκ)

    met.dgu[2, 2, 2] = (2*dᾱdκ*d2ᾱdκdᾱ*∇κ2 + dᾱdκ^2*d∇κ2dᾱ +
                        2*(d2ᾱdβdᾱ*dᾱdκ*∇κ∇β + dᾱdβ*d2ᾱdκdᾱ*∇κ∇β + dᾱdβ*dᾱdκ*d∇κ∇βdᾱ) +
                        2*dᾱdβ*∇β2*d2ᾱdβdᾱ + dᾱdβ^2*d∇β2dᾱ)

    met.dgu[2, 3, 1] = n0/R0^2 * ((d2ᾱdκ2*sqκ*sn*dn + dᾱdκ*dsqκ*sn*dn + dᾱdκ*sqκ*dsndκ*dn + dᾱdκ*sqκ*sn*ddndκ) +
                                  (d2ᾱdβdκ*dsqκ*dn*cn - 1/4*dᾱdβ*κ^(-3/2)*dn*cn + dᾱdβ*dsqκ*ddndκ*cn + dᾱdβ*dsqκ*dn*dcndκ))

    met.dgu[2, 3, 2] = n0/R0^2 * ((d2ᾱdκdᾱ*sqκ*sn*dn + dᾱdκ*sqκ*dsndᾱ*dn + dᾱdκ*sqκ*sn*ddndᾱ) +
                                  (d2ᾱdβdᾱ*dsqκ*dn*cn + dᾱdβ*dsqκ*ddndᾱ*cn + dᾱdβ*dsqκ*dn*dcndᾱ))


    met.dgu[3, 1, 1] = n0/R0^2 * (dsqκ*sn*dn + sqκ*dsndκ*dn + sqκ*sn*ddndκ)
    met.dgu[3, 1, 2] = n0/R0^2 * (sqκ * dsndᾱ*dn + sqκ*sn*ddndᾱ)
 
    met.dgu[3, 2, 1] = n0/R0^2 * ((d2ᾱdκ2*sqκ*sn*dn + dᾱdκ*dsqκ*sn*dn + dᾱdκ*sqκ*dsndκ*dn + dᾱdκ*sqκ*sn*ddndκ) +
                                  (d2ᾱdβdκ*dsqκ*dn*cn - 1/4*dᾱdβ*κ^(-3/2)*dn*cn + dᾱdβ*dsqκ*ddndκ*cn + dᾱdβ*dsqκ*dn*dcndκ))

    met.dgu[3, 2, 2] = n0/R0^2 * ((d2ᾱdκdᾱ*sqκ*sn*dn + dᾱdκ*sqκ*dsndᾱ*dn + dᾱdκ*sqκ*sn*ddndᾱ) +
                                  (d2ᾱdβdᾱ*dsqκ*dn*cn + dᾱdβ*dsqκ*ddndᾱ*cn + dᾱdβ*dsqκ*dn*dcndᾱ))

    met.gl .= inv(met.gu) 

    #using derivtive of inverse matrix
    #(K^-1)' = - K^-1 (K)' K^-1
    for i in 1:3
        met.dgl[:, :, i] =  -1 .* met.gl * met.dgu[:, :, i] * met.gl
    end
    #two ways of computing J are the same now!!! v nice.
    #met.J = sqrt(det(met.gl))# * R0 * sqrt(r2) #not certain about this!
    #same shape but different scale!!
    #this still does not match Axel's, off by factor of 2/m0, this could be explained by difference in α definition or κ??
    met.J[1] = K * w / (m0*π) * R0# *sqrt(r2))

    #J only function of κ
    met.dJ[1] = w * R0/(m0*π) * (E - (1-κ) * K) / (2*(1-κ)*κ)


end



function Axel_island_metric!(met::MetT, κ::Float64, ᾱ::Float64, φ::Float64, R0::Float64)

    #may need to change κ to be sqrt. Follow Axel's first.

    K, E = Elliptic.ellipke(κ)

    sn, cn, dn = Elliptic.ellipj(4.0*K / (2*π) * ᾱ, κ)

    F = Elliptic.F(4.0*K / (2*π) * ᾱ, κ)

    Z = Elliptic.E( 4.0*K / (2*π) * ᾱ, κ) - E/K * F


    #hopefully we don't actually need the island for this...
    #think we will, probably need an island submodule... big RIP for combining everything.
    #maybe we need an island coords flag or something to separate them all, or just an island grid, which is the same functionally,
    #but uses κ etc and has a metric function that accepts an island.
    r02 = (0.5)^2
    w = 0.05
    m0 = 2
    n0 = -1
    r2 = r02 + w * κ * cn

    #β is used the same as Axel, but we have replaced βs with ᾱ to match our eventual method.

    #this assumes cylindrical metric elements. Hence the r^2 and R0^2, for torus they are more complicated.
    ∇κ2 = r02/w^2 * (4.0*cn^2 * r2/r02 + w^2/r02 * 0.25 * sn^2 * dn^2 * m0^2/r2 + n0^2/R0^2)

    ∇β2 = r02 / (w^2*κ^2) * (4.0 * sn^2 * r2/r02 + w^2/r02 * 0.25 * cn^2*dn^2 * m0^2/r2 + n0^2/R0^2)

    ∇κ∇β = r02 / (w^2 * κ) * (-4.0 * sn*cn * r2/r02 + w^2/r02 * 0.25 * cn^2*dn^2 * m0^2/r2 + n0^2/R0^2)
    
    dᾱdβ = 2 * π / (4*K) / dn

    #this term does not match paper, lone κ should be on top.
    dᾱdκ = 2 * π / (4*K * κ * (1-κ^2)) * (Z - κ^2*cn*sn / dn)

    met.J[1] = 4*R0 * κ * K * w / (2*π*m0^2)
    
    
    met.gu[1, 1] = ∇κ2
    met.gu[1, 2] = dᾱdβ * ∇κ∇β +  dᾱdκ * ∇κ2
    met.gu[1, 3] = 0.5 * n0/R0^2 * sn*dn

    met.gu[2, 1] = dᾱdβ * ∇κ∇β +  dᾱdκ * ∇κ2
    met.gu[2, 2] = dᾱdβ^2 * ∇β2 + 2 * dᾱdβ * dᾱdκ * ∇κ∇β + dᾱdκ^2 * ∇κ2
    met.gu[2, 3] = dᾱdβ * 0.5*n0/R0^2*cn*dn/κ + dᾱdκ*0.5*n0/R0^2*sn*dn

    met.gu[3, 1] = 0.5 * n0/R0^2 * sn*dn
    met.gu[3, 2] = dᾱdβ * 0.5*n0/R0^2*cn*dn/κ + dᾱdκ*0.5*n0/R0^2*sn*dn
    met.gu[3, 3] = 1/R0^2


    met.gl = inv(met.gu) 


end


"""
    slab_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)

Function that fills out the MetT structure for a slab of plasma/cartesian geometry.
"""
function slab_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)

    met.gl[1, 1] = 1.0
    met.gl[2, 2] = 1.0
    met.gl[3, 3] = R0^2

    met.gu[1, 1] = 1.0
    met.gu[2, 2] = 1.0
    met.gu[3, 3] = 1 / R0^2

    met.J[1] = R0


end



"""
    no_delta_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)

Function that fills out the MetT struct with Δ'=0, otherwise identical to toroidal_metric. Used for comparison with literature.
"""
function no_delta_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)
    
    
    Δp = 0.0
    Δpp = 0.0

    ϵ = r/R0

    η = 1/2 * (ϵ+Δp)

    ηp = 1/2* (1/R0 + Δpp)

    met.J[1] = r * R0 * (1+2*ϵ*cos(θ))

    met.gl[1, 1] = 1-2*Δp * cos(θ)
    met.gl[1, 2] = r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gl[2, 1] = r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gl[2, 2] = r^2*(1+4*η*cos(θ) + 4*η^2)
    met.gl[3, 3] = R0^2*(1+2*ϵ*cos(θ))


    met.gu[1, 1] = 1+2*Δp * cos(θ)
    met.gu[1, 2] = -1/r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gu[2, 1] = -1/r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gu[2, 2] = 1/r^2*(1-2*(ϵ+Δp)*cos(θ))
    met.gu[3, 3] = 1/R0^2*(1-2*ϵ*cos(θ))

    
    met.dJ[1] = R0 + 4*r * cos(θ)
    met.dJ[2] = -2 * r * R0*ϵ * sin(θ)

    #first two indicies give metric element, while third is derivative,
    #eg [1, 2, 3] is ∂g_{12}/∂ζ
    met.dgl[1, 1, 1] = -2*Δpp * cos(θ)
    met.dgl[1, 1, 2] = 2*Δp * sin(θ)

    met.dgl[1, 2, 1] = ((ϵ + Δp + r*Δpp) + r*(1/R0 + 2*Δpp)) * sin(θ)
    met.dgl[1, 2, 2] = r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgl[2, 1, 1] = ((ϵ + Δp + r*Δpp) + r*(1/R0 + 2*Δpp)) * sin(θ)
    met.dgl[2, 1, 2] = r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgl[2, 2, 1] = 2*r*(1+4*η*cos(θ) + 4*η^2) + r^2 * (4*ηp*cos(θ) + 8*η * ηp)
    met.dgl[2, 2, 2] = -r^2*(4*η*sin(θ))

    met.dgl[3, 3, 1] = 2*R0*cos(θ)
    met.dgl[3, 3, 2] = -2*R0^2*ϵ*sin(θ)

    met.dgu[1, 1, 1] = 2*Δpp * cos(θ)
    met.dgu[1, 1, 2] = -2*Δp * sin(θ)

    met.dgu[1, 2, 1] = (1/r^2 * (ϵ + Δp + r*Δpp) - 1/r*(1/R0 + 2*Δpp)) * sin(θ)
    met.dgu[1, 2, 2] = -1/r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgu[2, 1, 1] = (1/r^2 * (ϵ + Δp + r*Δpp) - 1/r*(1/R0 + 2*Δpp)) * sin(θ)
    met.dgu[2, 1, 2] = -1/r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgu[2, 2, 1] = -2/r^3 * (1-2*(ϵ+Δp)*cos(θ)) + 1/r^2 * (-2*(1/R0+Δpp)*cos(θ))
    met.dgu[2, 2, 2] = 2/r^2*(ϵ+Δp)*sin(θ)

    met.dgu[3, 3, 1] = -2*cos(θ)/R0^3
    met.dgu[3, 3, 2] = 2*ϵ*sin(θ)/R0^2


end

"""
    diagonal_toroidal_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)

Toroidal metric retaining only the diagonal elements. Used for comparison with two mode models.
"""
function diagonal_toroidal_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)
    #TODO

    #this is cooked in its current form.
    #gives a completly different frequency.
    #pretty cooked, :(
    
    Δp = r/(4*R0)
    Δpp = 1/(4*R0)

    ϵ = r/R0

    η = 1/2 * (ϵ+Δp)

    ηp = 1/2* (1/R0 + Δpp)

    met.J[1] = r * R0 * (1+2*ϵ*cos(θ))

    #don't really know anyting about gl from paper.
    met.gl[1, 1] = 1-2*Δp * cos(θ)
    #met.gl[1, 2] = r*(ϵ + Δp + r*Δpp) * sin(θ)
    #met.gl[2, 1] = r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gl[2, 2] = r^2*(1+4*η*cos(θ) + 4*η^2)
    met.gl[3, 3] = R0^2*(1+2*ϵ*cos(θ)) #this term seems to be difference between 77 vs 12. hmmm.


    met.gu[1, 1] = 1+2*Δp * cos(θ)
    #met.gu[1, 2] = -1/r*(ϵ + Δp + r*Δpp) * sin(θ)
    #met.gu[2, 1] = -1/r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gu[2, 2] = 1/r^2*(1-2*(ϵ+Δp)*cos(θ))
    met.gu[3, 3] = 1/R0^2*(1-2*ϵ*cos(θ)) 

    
    met.dJ[1] = R0 + 4*r * cos(θ)
    #met.dJ[2] = -2 * r * R0*ϵ * sin(θ)

    #first two indicies give metric element, while third is derivative,
    #eg [1, 2, 3] is ∂g_{12}/∂ζ
    met.dgl[1, 1, 1] = -2*Δpp * cos(θ)
    #met.dgl[1, 1, 2] = 2*Δp * sin(θ)

    #met.dgl[1, 2, 1] = ((ϵ + Δp + r*Δpp) + r*(1/R0 + 2*Δpp)) * sin(θ)
    #met.dgl[1, 2, 2] = r*(ϵ + Δp + r*Δpp) * cos(θ)

    #met.dgl[2, 1, 1] = ((ϵ + Δp + r*Δpp) + r*(1/R0 + 2*Δpp)) * sin(θ)
    #met.dgl[2, 1, 2] = r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgl[2, 2, 1] = 2*r#*(1+4*η*cos(θ) + 4*η^2) + r^2 * (4*ηp*cos(θ) + 8*η * ηp)
    #met.dgl[2, 2, 2] = -r^2*(4*η*sin(θ))

    #met.dgl[3, 3, 1] = 2*R0*cos(θ)
    #met.dgl[3, 3, 2] = -2*R0^2*ϵ*sin(θ)

    met.dgu[1, 1, 1] = 2*Δpp * cos(θ)
    #met.dgu[1, 1, 2] = -2*Δp * sin(θ)

    #met.dgu[1, 2, 1] = (1/r^2 * (ϵ + Δp + r*Δpp) - 1/r*(1/R0 + 2*Δpp)) * sin(θ)
    #met.dgu[1, 2, 2] = -1/r*(ϵ + Δp + r*Δpp) * cos(θ)

    #met.dgu[2, 1, 1] = (1/r^2 * (ϵ + Δp + r*Δpp) - 1/r*(1/R0 + 2*Δpp)) * sin(θ)
    #met.dgu[2, 1, 2] = -1/r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgu[2, 2, 1] = -2/r^3 #* (1-2*(ϵ+Δp)*cos(θ)) + 1/r^2 * (-2*(1/R0+Δpp)*cos(θ))
    #met.dgu[2, 2, 2] = 2/r^2*(ϵ+Δp)*sin(θ)

    #met.dgu[3, 3, 1] = -2*cos(θ)/R0^3
    #met.dgu[3, 3, 2] = 2*ϵ*sin(θ)/R0^2


end



"""
    cylindrical_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)

Cylindrical limit of toroidal metric, equivalent to taking R0→∞.
"""
function rad_cylindrical_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)
    #this is regular old cylindrical for comparing our weak form

    met.J[1] = r * R0

    met.gl[1, 1] = 1

    met.gl[2, 2] = r^2

    met.gl[3, 3] = R0^2


    met.gu[1, 1] = 1

    met.gu[2, 2] = 1/r^2

    met.gu[3, 3] = 1/R0^2


    met.dgl[2, 2, 1] = 2*r

    met.dJ[1] = R0
    
    met.dgu[2, 2, 1] = -2 / r^3

end

function flux_cylindrical_metric!(met::MetT, ψ::Float64, θ::Float64, ζ::Float64, R0::Float64)
        #this is regular old cylindrical for comparing our weak form

    r = sqrt(2*ψ)
    dψdr = r 

    met.J[1] = R0

    met.gl[1, 1] = 1/r^2

    met.gl[2, 2] = r^2

    met.gl[3, 3] = R0^2


    met.gu[1, 1] = r^2

    met.gu[2, 2] = 1/r^2

    met.gu[3, 3] = 1/R0^2


    met.dgl[1, 1, 1] = -2 / r^3 / dψdr
    met.dgl[2, 2, 1] = 2*r / dψdr

    met.dgu[1, 1, 1] = 2*r / dψdr
    met.dgu[2, 2, 1] = -2 / r^3 / dψdr

end



"""
    flux_toroidal_metric!(met::MetT, ψ::Float64, θ::Float64, ζ::Float64, R0::Float64)

Function for toroidal metric with flux as the radial coordinate. Used by island continuum. 
Currently only computes only what is required for island continuum.
"""
#no idea if this is actually being used, if so it should be in MIDIslands tbh
function isl_flux_toroidal_metric!(met::MetT, ψ::Float64, θ::Float64, ζ::Float64, R0::Float64)


    r = sqrt(2*ψ) #B0=1
    dψdr = r 

    d2ψdr2 = 1 #simplest way to use previous results!

    Δp = r/(4*R0)
    Δpp = 1/(4*R0)

    ϵ = r/R0

    η = 1/2 * (ϵ+Δp)

    ηp = 1/2* (1/R0 + Δpp)

    met.J[1] = r * R0 * (1+2*ϵ*cos(θ)) / dψdr

    #guessing that gl would be divided by dψdr, will have to compare the resulting metrics.
    met.gl[1, 1] = (1-2*Δp * cos(θ)) / dψdr^2
    met.gl[1, 2] = r*(ϵ + Δp + r*Δpp) * sin(θ) / dψdr
    met.gl[2, 1] = r*(ϵ + Δp + r*Δpp) * sin(θ) / dψdr
    met.gl[2, 2] = r^2*(1+4*η*cos(θ) + 4*η^2)
    met.gl[3, 3] = R0^2*(1+2*ϵ*cos(θ))


    met.gu[1, 1] = (1+2*Δp * cos(θ)) * dψdr^2
    met.gu[1, 2] = -1/r*(ϵ + Δp + r*Δpp) * sin(θ) * dψdr
    met.gu[2, 1] = -1/r*(ϵ + Δp + r*Δpp) * sin(θ) * dψdr
    met.gu[2, 2] = 1/r^2*(1-2*(ϵ+Δp)*cos(θ))
    met.gu[3, 3] = 1/R0^2*(1-2*ϵ*cos(θ))

end




"""
    flux_toroidal_metric!(met::MetT, ψ::Float64, θ::Float64, ζ::Float64, R0::Float64)

Function for toroidal metric with flux as the radial coordinate.
"""
function flux_toroidal_metric!(met::MetT, ψ::Float64, θ::Float64, ζ::Float64, R0::Float64)


    #actually tested this now, and looks good.

    #this will be assuming B0=1 everywhere. Not sure if it ever needs to be changed.
    r = sqrt(2*ψ) #B0=1
    dψdr = r 

    #drdψ = 1 / r

    #d2ψdr2 = 1 #simplest way to use previous results!

    Δp = r/(4*R0)
    Δpp = 1/(4*R0)

    ϵ = r/R0

    η = 1/2 * (ϵ+Δp)

    ηp = 1/2* (1/R0 + Δpp)

    #so we will be using the radial form and modifying for flux.

    met.J[1] = R0 * (1+2*ϵ*cos(θ))

    #guessing that gl would be divided by dψdr, will have to compare the resulting metrics.
    #met.gl[1, 1] = (1-2*Δp * cos(θ)) / dψdr^2
    #met.gl[1, 2] = r*(ϵ + Δp + r*Δpp) * sin(θ) / dψdr
    #met.gl[2, 1] = r*(ϵ + Δp + r*Δpp) * sin(θ) / dψdr
    #met.gl[2, 2] = r^2*(1+4*η*cos(θ) + 4*η^2)
    #met.gl[3, 3] = R0^2*(1+2*ϵ*cos(θ))


    met.gu[1, 1] = r^2 * (1+2*Δp * cos(θ)) 
    #met.gu[1, 2] = -1/r*(ϵ + Δp + r*Δpp) * sin(θ) * dψdr
    #met.gu[2, 1] = -1/r*(ϵ + Δp + r*Δpp) * sin(θ) * dψdr
    met.gu[1, 2] = -(ϵ + Δp + r*Δpp) * sin(θ)
    met.gu[2, 1] = -(ϵ + Δp + r*Δpp) * sin(θ) 
    met.gu[2, 2] = 1/r^2*(1-2*(ϵ+Δp)*cos(θ))
    met.gu[3, 3] = 1/R0^2*(1-2*ϵ*cos(θ))

    met.gl .= inv(met.gu)


    met.dJ[1] = 2*cos(θ) / dψdr
    met.dJ[2] = -R0 * 2*ϵ*sin(θ)


    met.dgu[1, 1, 1] = (2*r * (1+2*Δp * cos(θ)) + r^2*(2*Δpp*cos(θ))) / dψdr
    met.dgu[1, 1, 2] = r^2 * (-2*Δp * sin(θ)) 

    met.dgu[1, 2, 1] = -(1/R0 + 2*Δpp) * sin(θ) / dψdr
    met.dgu[1, 2, 2] = -(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgu[2, 1, 1] = -(1/R0 + 2*Δpp) * sin(θ) / dψdr
    met.dgu[2, 1, 2] = -(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgu[2, 2, 1] = (-2/r^3 * (1-2*(ϵ+Δp)*cos(θ)) + 1/r^2*(-2*(1/R0 + Δpp)*cos(θ))) / dψdr
    met.dgu[2, 2, 2] = -1/r^2*(-2*(ϵ+Δp)*sin(θ)) 

    met.dgu[3, 3, 1] = 1/R0^2*(-2/R0*cos(θ)) / dψdr
    met.dgu[3, 3, 2] = -1/R0^2*(-2*ϵ*sin(θ))

    for i in 1:3
        met.dgl[:, :, i] = - met.gl * met.dgu[:, :, i] * met.gl
    end

end
