"""

Module for numerical integration using Gaussian quadrature.
"""
module Integration


import FastGaussQuadrature: gausslegendre


using ..Grids


include("GaussQuadrature.jl")

export gauss_integrate, gauss_points


end
