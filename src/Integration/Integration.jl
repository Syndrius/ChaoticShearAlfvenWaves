"""

Module for numerical integration using Gaussian quadrature.
"""
module Integration

#perhaps this should be called Gauss quadrature not integration!
using FastGaussQuadrature
using ..Grids


include("GaussQuadrature.jl")

export gauss_integrate, gauss_points


end
