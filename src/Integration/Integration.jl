"""

Module for numerical integration using Gaussian quadrature.
"""
module Integration


using FastGaussQuadrature

using ..Grids


#good
include("GaussQuadrature.jl")

export gauss_integrate, gauss_points


end
