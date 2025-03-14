
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
using MID.QFM #not ideal, may need to move the QFM structure into Structures.

using Plots
using LaTeXStrings
using Printf
using FFTW


include("PotentialPlot.jl")

export potential_plot


include("ContourPlot.jl")

export contour_plot
export contour_zeta_plot


include("SurfacePlot.jl")

export surface_plot


include("ContinuumPlot.jl")

export continuum_plot


include("QFMPlotting.jl")

#name of this is unfort, may need to change to qfm_surf_plot or something.
export plot_surfs

end