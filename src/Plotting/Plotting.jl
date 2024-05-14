
"""
May want to change this to post-processing or something!

#this module still hasn't been fixed up yet!
#this is in a weird spot as much of the functionality is in MIDViz.
#not sure exactly what we want in here!
#this is still ~90% of the compilation, could be worth removing all plotting from MID

"""

module Plotting

using Plots
#Plots.set_default_backend!(:plotlyjs)
#we probbaly want to use the plotly backend for inspection of specific frequencies.
#using PlotlyJS

include("PotentialPlot.jl")

export plot_potential
export find_ind

include("ContinuumPlot.jl")

export plot_continuum
export reconstruct_continuum


end