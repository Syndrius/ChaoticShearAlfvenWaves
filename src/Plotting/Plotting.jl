
"""
May want to change this to post-processing or something!

#this module still hasn't been fixed up yet!
#this is in a weird spot as much of the functionality is in MIDViz.
#not sure exactly what we want in here!
#this is still ~90% of the compilation, could be worth removing all plotting from MID

"""

module Plotting

using MID.Inputs
#using MID.

using Plots
using LaTeXStrings
using Printf
using FFTW
#Plots.set_default_backend!(:plotlyjs)
#we probbaly want to use the plotly backend for inspection of specific frequencies.
#using PlotlyJS

include("PotentialPlot.jl")

export plot_potential
export plot_sum_potential
export find_ind
export plot_phi_surface
export construct_surface

include("ContinuumPlot.jl")

export plot_continuum
export reconstruct_continuum
export reconstruct_continuum_n
export mode_structure


end