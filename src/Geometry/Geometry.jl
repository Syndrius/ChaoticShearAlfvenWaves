"""

Module for the computing the metric tensor in different coordinate systems.
"""
module Geometry


using ..Structures


using Elliptic
import LinearAlgebra: inv


export init_geometry


include("Toroidal.jl")

export radial_toroidal_metric!, flux_toroidal_metric!


include("Cylindrical.jl")

export radial_cylindrical_metric!, flux_cylindrical_metric!


include("Island.jl")

export island_metric!


"""
    init_geometry(type::Symbol=:tor; R0::Float64=4.0)

Initialises the geometry struct. This struct is recreated when the ProblemT struct is created as the geometry can depend on the fields used.
"""
function init_geometry(type::Symbol=:tor; R0::Float64=4.0)

    #flux_toroidal_metric! is a placeholder.
    return GeometryT(type, R0, flux_toroidal_metric!, 1.0, 1.0)
end


end
