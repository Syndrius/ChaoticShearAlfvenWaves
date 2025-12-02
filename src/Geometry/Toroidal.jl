"""
    flux_toroidal_metric!(met::MetT, ψ::Float64, θ::Float64, φ::Float64, R0::Float64)

Function for toroidal metric with flux as the radial coordinate.
"""
function flux_toroidal_metric!(met::MetT, ψ::Float64, θ::Float64, φ::Float64, R0::Float64)


    r = sqrt(2*ψ) #B0=1
    dψdr = r 

    Δp = r/(4*R0)
    Δpp = 1/(4*R0)

    ϵ = r/R0

    η = 1/2 * (ϵ+Δp)

    ηp = 1/2* (1/R0 + Δpp)

    met.J[1] = R0 * (1+2*ϵ*cos(θ))


    met.gu[1, 1] = r^2 * (1+2*Δp * cos(θ)) 
    met.gu[1, 2] = -(ϵ + Δp + r*Δpp) * sin(θ)
    met.gu[2, 1] = -(ϵ + Δp + r*Δpp) * sin(θ) 
    met.gu[2, 2] = 1/r^2*(1-2*(ϵ+Δp)*cos(θ))
    met.gu[3, 3] = 1/R0^2*(1-2*ϵ*cos(θ))

    met.gl .= inv(met.gu)

    met.dJ[1] = 2*cos(θ) / dψdr
    met.dJ[2] = -R0 * 2*ϵ*sin(θ)


    #first two indicies give metric element, while third is derivative,
    #eg [1, 2, 3] is ∂g_{12}/∂x3
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


"""
    radial_toroidal_metric!(met::MetT, r::Float64, θ::Float64, φ::Float64, R0::Float64)

Function that fills out the MetT struct for toroidal geometry. Metric elements taken from Energetic Particles in Tokamak Plasmas by Sergai Sharapov. Straight field line coordinates are radius (r), generalised poloidal angle (θ) and generalised toroidal angle (φ), equal to negative of true toroidal angle. Additionally we assume low shear and approximate Δ' ≈ r/(4*R0).
"""
function radial_toroidal_metric!(met::MetT, r::Float64, θ::Float64, φ::Float64, R0::Float64)
    

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
    #eg [1, 2, 3] is ∂g_{12}/∂φ
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

