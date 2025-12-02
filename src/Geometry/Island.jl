"""
    island_metric!(met::MetT, κ::Float64, ᾱ::Float64, φ::Float64, R0::Float64, isl::CoordIslandT)

Function that fills out the MetT structure for magnetic island striaght field line geometry.
Derived follwing Konies et al 2024, using the notation of Qu and Hole 2023.
κ is the radial variable, measure from island center, ᾱ is the striaghtened helical angle and τ is a typical toroidal angle.
This assumed the original geometry was cylindrical.
"""
function island_metric!(met::MetT, κ::Float64, ᾱ::Float64, τ::Float64, R0::Float64, isl::CoordIslandT)

    #first we compute our intermediate metric (κ, β, τ)
    #where √κsin(β) = sin(m0 α/2)

    #we have, for the first transformation
    #κ = 4/w^2 * (ψ-ψ0)^2 + sin^2(m0*α/2) -> κ(ψ, θ, φ)
    #β = arctan(w/2 * sin(m0*α/2)/(ψ-ψ0)) -> β(ψ, θ, φ)
    #τ = φ -> τ(ψ, θ, φ)

    #We then straighten these coordinates via
    #ᾱ = π/(2*K(κ)) * F(β, κ) -> ᾱ(κ, β, τ)

    #other relations
    #ψ-ψ0 = w*sqrt(κ)/2 * cos(β)
    #but we will use β as an intermediate for derivatives to reduce complexity.

    K, E = ellipke(κ)

    #shorhand for the argument in the coordinate transformation
    arg = 2*K / π * ᾱ

    sn, cn, dn = ellipj(arg, κ)
    β = Jacobi.am(arg, κ)

    Eβ = Elliptic.E(β, κ) #equiv to JacobiEpsilon
    Fβ = Elliptic.F(β, κ)

    Z = Eβ - E / K * Fβ

    #extract the island info
    ψ0 = isl.ψ0
    w = isl.w
    m0 = isl.m0
    n0 = isl.n0

    ψ = ψ0 + w * sqrt(κ) / 2  * cn

    sqκ = sqrt(κ)
    dsqκ = 1/(2*sqκ)

    chat = dn 

    bhat = (m0^2/(2ψ) + n0^2 / R0^2) * chat^2

    ∇κ∇β = (-16*ψ/w^2 + bhat/2) * cn * sn


    ∇β2 = 8*ψ*sn^2 / (w^2*κ) + bhat/(4κ) * cn^2

    ∇β∇τ = n0/(2*sqrt(κ)*R0^2) * cn * dn #/ R0^2 #this term is probbaly negligible.

    ∇κ2 = 32 * ψ * κ / w^2 * cn^2 + bhat * κ * sn^2

    dᾱdκ = π / (4*κ * (1-κ) * K) * (Z - κ*sn*cn / dn)

    dᾱdβ = π / (2 * K * chat)

    ∇κ∇ᾱ = dᾱdκ * ∇κ2 + dᾱdβ * ∇κ∇β

    ∇κ∇τ = n0/R0^2 * sqrt(κ) * sn * dn

    ∇ᾱ2 = (dᾱdκ)^2 * ∇κ2 + 2 * dᾱdκ * dᾱdβ * ∇κ∇β + (dᾱdβ)^2 * ∇β2

    ∇ᾱ∇τ = dᾱdκ * ∇κ∇τ + dᾱdβ * ∇β∇τ

    ∇τ2 = 1 / R0^2

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

    met.J[1] = w/(m0*π) * K * R0

    #now the derivatives 

    cd = cn / dn
    sc = sn / cn

    dKdκ = (E - (1-κ)*K) / (2*(1-κ)*κ)
    dEdκ = (E - K)/(2*κ)
    dargdκ = 2/π * ᾱ * dKdκ
    dargdᾱ = 2/π * K
    dβdκ = 1/(2*(κ-1)*κ) * (dn * ((κ-1)*arg + Eβ) - κ*cn * sn) + dn * dargdκ
    dβdᾱ = dn * dargdᾱ

    dEβdκ = (Eβ - Fβ)/(2*κ) + sqrt(1-κ*sin(β)^2)*dβdκ
    dEβdᾱ = sqrt(1-κ*sin(β)^2)*dβdᾱ
    dFβdκ = 1/4 * (2*Eβ/(κ-κ^2) - 2*Fβ/κ + sin(2*β) / ((κ-1)*sqrt(1-κ*sin(β)^2)) + 4*dβdκ/sqrt(1-κ*sin(β)^2))
    dFβdᾱ = dβdᾱ / sqrt(1-κ*sin(β)^2)

    dZdκ = dEβdκ - (dEdκ / K * Fβ - E/K^2 * dKdκ * Fβ + E/K * dFβdκ)
    dZdᾱ = dEβdᾱ - (E/K * dFβdᾱ)


    ddndκ = 1/2*cn * sn * (1/(1-κ)*((κ-1) * arg + Eβ - dn * sc) - 2* κ*dargdκ)
    ddndᾱ = -κ * cn * sn * dargdᾱ


    dsndκ = 1/(2*(κ-1)*κ) * cn * dn * ((κ-1)*arg + Eβ - κ * cd*sn + 2*(κ-1) * κ*dargdκ)
    dsndᾱ = cn * dn * dargdᾱ

    dcndκ = 1/(2*(κ-1)*κ) * dn * sn * (-((κ-1)*arg) - Eβ + κ * cd*sn - 2*(κ-1) * κ*dargdκ)
    dcndᾱ = -dn * sn * dargdᾱ

    dψdκ = w * 1/2*(κ)^(-1/2)/2*cn + w*sqrt(κ)/2 * dcndκ
    dψdᾱ = w * sqκ / 2 * dcndᾱ

    dbhatdκ = -m0^2/(2ψ^2) * dψdκ * dn^2 + (m0^2/(2ψ) + n0^2 / R0^2) * 2*dn*ddndκ
    dbhatdᾱ = -m0^2/(2ψ^2) * dψdᾱ * dn^2 + (m0^2/(2ψ) + n0^2 / R0^2) * 2*dn*ddndᾱ

    d2ᾱdκ2 = (-π / (4*κ*(1-κ)*K)^2 * 4*((1-κ)*K - κ*K + κ*(1-κ)*dKdκ) * (Z-κ*sn*cn / dn)
              + π / (4*κ*(1-κ)*K) * (dZdκ - (sn*cn/dn + κ*dsndκ*cn/dn + κ*sn*dcndκ/dn - κ*sn*cn/dn^2*ddndκ)))

    d2ᾱdβdκ = -π / (2*K*dn)^2 * (2*dKdκ*dn + 2*K*ddndκ)

    d2ᾱdβdᾱ = -π / (2 * K * dn)^2 * (2*K * ddndᾱ)

    d2ᾱdκdᾱ = π / (4*κ * (1-κ) * K) * (dZdᾱ - κ*(dsndᾱ*cn/dn + sn*dcndᾱ/dn - sn*cn/dn^2 * ddndᾱ))

    d∇κ∇βdκ = (-16*dψdκ/w^2 + dbhatdκ/2) * cn * sn + (-16*ψ/w^2 + bhat/2) * (dcndκ*sn + cn * dsndκ)

    d∇κ∇βdᾱ = (-16 * dψdᾱ / w^2 + dbhatdᾱ/2) * cn * sn + (-16*ψ/w^2 + bhat/2) * (dcndᾱ*sn + cn * dsndᾱ)

    d∇β2dκ = (8/w^2*(dψdκ * sn^2 / κ + ψ * 2 * sn * dsndκ / κ - ψ*sn^2 / κ^2) 
              + dbhatdκ /(4*κ) * cn^2 - bhat / (4κ^2) * cn^2 + bhat /(4*κ) * 2*cn*dcndκ)

    d∇β2dᾱ = 8/(w^2*κ)*(dψdᾱ * sn^2 + ψ * 2 * sn * dsndᾱ) + dbhatdᾱ /(4*κ) * cn^2 + bhat /(4*κ) * 2*cn*dcndᾱ

    d∇β∇τdκ = n0 / (2*R0^2) * (-dsqκ/κ * cn * dn + dcndκ*dn/sqκ + cn*ddndκ/sqκ)
    d∇β∇τdᾱ = n0 / (2*R0^2) * (dcndᾱ*dn/sqκ + cn*ddndᾱ/sqκ)

    
    met.dJ[1] = R0 * w / (2*m0*π) * (E - (1-κ) * K) / ((1-κ)*κ)

    d∇κ2dκ = 32 / w^2 * (dψdκ*κ*cn^2 + ψ*cn^2 + ψ*κ*2*cn*dcndκ) + (dbhatdκ*κ*sn^2 + bhat*sn^2 + bhat*κ*2*sn*dsndκ)
    d∇κ2dᾱ = 32 * κ / w^2 * (dψdᾱ * cn^2 + ψ*2*cn*dcndᾱ) + κ*(bhat*2*sn*dsndᾱ + dbhatdᾱ *sn^2)

    d∇κ∇ᾱdκ = d2ᾱdκ2 * ∇κ2 + dᾱdκ * d∇κ2dκ + d2ᾱdβdκ * ∇κ∇β + dᾱdβ * d∇κ∇βdκ
    d∇κ∇ᾱdᾱ = d2ᾱdκdᾱ * ∇κ2 + dᾱdκ * d∇κ2dᾱ + d2ᾱdβdᾱ * ∇κ∇β + dᾱdβ * d∇κ∇βdᾱ

    d∇κ∇τdκ = n0/R0^2 * (dsqκ * sn * dn + sqκ * dsndκ * dn + sqκ * sn * ddndκ)
    d∇κ∇τdᾱ = n0/R0^2 * (sqκ * dsndᾱ * dn + sqκ * sn * ddndᾱ)

    d∇ᾱ2dκ = (2*(dᾱdκ)*d2ᾱdκ2 * ∇κ2 + dᾱdκ^2 * d∇κ2dκ 
              + 2 * (d2ᾱdκ2 * dᾱdβ * ∇κ∇β + dᾱdκ * d2ᾱdβdκ * ∇κ∇β + dᾱdκ * dᾱdβ * d∇κ∇βdκ) 
              + 2*(dᾱdβ)*d2ᾱdβdκ * ∇β2 + dᾱdβ^2 * d∇β2dκ)

    d∇ᾱ2dᾱ = (2*(dᾱdκ)*d2ᾱdκdᾱ * ∇κ2 + dᾱdκ^2 * d∇κ2dᾱ
              + 2 * (d2ᾱdκdᾱ * dᾱdβ * ∇κ∇β + dᾱdκ * d2ᾱdβdᾱ * ∇κ∇β + dᾱdκ * dᾱdβ * d∇κ∇βdᾱ) 
              + 2*(dᾱdβ)*d2ᾱdβdᾱ * ∇β2 + dᾱdβ^2 * d∇β2dᾱ)

    d∇ᾱ∇τdκ = d2ᾱdκ2 * ∇κ∇τ + dᾱdκ * d∇κ∇τdκ + d2ᾱdβdκ * ∇β∇τ + dᾱdβ * d∇β∇τdκ
    d∇ᾱ∇τdᾱ = d2ᾱdκdᾱ * ∇κ∇τ + dᾱdκ * d∇κ∇τdᾱ + d2ᾱdβdᾱ * ∇β∇τ + dᾱdβ * d∇β∇τdᾱ


    met.dgu[1, 1, 1] = d∇κ2dκ
    met.dgu[1, 1, 2] = d∇κ2dᾱ
    met.dgu[1, 1, 3] = 0.0

    met.dgu[1, 2, 1] = d∇κ∇ᾱdκ 
    met.dgu[1, 2, 2] = d∇κ∇ᾱdᾱ 
    met.dgu[1, 2, 3] = 0.0 

    met.dgu[1, 3, 1] = d∇κ∇τdκ 
    met.dgu[1, 3, 2] = d∇κ∇τdᾱ 
    met.dgu[1, 3, 3] = 0.0 

    met.dgu[2, 1, 1] = d∇κ∇ᾱdκ 
    met.dgu[2, 1, 2] = d∇κ∇ᾱdᾱ 
    met.dgu[2, 1, 3] = 0.0

    met.dgu[2, 2, 1] = d∇ᾱ2dκ 
    met.dgu[2, 2, 2] = d∇ᾱ2dᾱ 
    met.dgu[2, 2, 3] = 0.0

    met.dgu[2, 3, 1] = d∇ᾱ∇τdκ 
    met.dgu[2, 3, 2] = d∇ᾱ∇τdᾱ 
    met.dgu[2, 3, 3] = 0.0

    met.dgu[3, 1, 1] = d∇κ∇τdκ 
    met.dgu[3, 1, 2] = d∇κ∇τdᾱ 
    met.dgu[3, 1, 3] = 0.0

    met.dgu[3, 2, 1] = d∇ᾱ∇τdκ 
    met.dgu[3, 2, 2] = d∇ᾱ∇τdᾱ 
    met.dgu[3, 2, 3] = 0.0

    met.dgu[3, 3, 1] = 0.0
    met.dgu[3, 3, 2] = 0.0
    met.dgu[3, 3, 3] = 0.0 

    #using derivtive of inverse matrix
    #(K^-1)' = - K^-1 (K)' K^-1
    for i in 1:3
        met.dgl[:, :, i] =  -1 .* met.gl * met.dgu[:, :, i] * met.gl
    end

end 
