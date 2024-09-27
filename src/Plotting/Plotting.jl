
"""

Contains functions for plotting the output potential and frequencies, includes
 - Plotting the shear Alfven continuum.
 - Plotting the mode structure as a function of radius.
 - Plotting r, θ contours of the potential.
 - Plotting the r, θ surface of the potential.
 
Additional plotting functions can be found in the separate package MIDViz.jl.
"""

module Plotting

using MID.Structures

using Plots
using LaTeXStrings
using Printf
using FFTW


include("PotentialPlot.jl")

export potential_plot


include("ContourPlot.jl")

export contour_plot


include("SurfacePlot.jl")

export surface_plot


include("ContinuumPlot.jl")

export continuum_plot


end