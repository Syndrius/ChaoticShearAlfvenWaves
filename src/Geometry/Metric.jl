#this file should be split amongst different files for each type of geometry.
#may still want to keep this main file contianing the struct, although we may eventually want to move that to Structures.

#=
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
=#

#template metric function
function metric!(met::MetT, x1::Float64, x2::Float64, x3::Float64, R0::Float64)
end

"""
    toroidal_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)

Function that fills out the MetT struct for toroidal geometry. Metric elements taken from Energetic Particles in Tokamak Plasmas by Sergai Sharapov. Straight field line coordinates are radius (r), generalised poloidal angle (θ) and generalised toroidal angle (ζ), equal to negative of true toroidal angle. Additionally we assume low shear and approximate Δ' ≈ r/(4*R0).
"""
function radial_toroidal_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)
    
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
function radial_cylindrical_metric!(met::MetT, r::Float64, θ::Float64, ζ::Float64, R0::Float64)
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
