"""
    flux_cylindrical_metric!(met::MetT, ψ::Float64, θ::Float64, φ::Float64, R0::Float64)

Computes the metric tensor in the cylindrical limit with toroidal flux as the radial coordinate.
"""
function flux_cylindrical_metric!(met::MetT, ψ::Float64, θ::Float64, φ::Float64, R0::Float64)

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
    radial_cylindrical_metric!(met::MetT, r::Float64, θ::Float64, φ::Float64, R0::Float64)

Cylindrical limit of toroidal metric, equivalent to taking R0→∞.
"""
function radial_cylindrical_metric!(met::MetT, r::Float64, θ::Float64, φ::Float64, R0::Float64)

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

