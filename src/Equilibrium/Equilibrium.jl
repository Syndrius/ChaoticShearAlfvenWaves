"""

Module for the equilibrium state of the plasma. This includes the density and the magnetic field, computed from an input q-profile and input magnetic islands.
"""
module Equilibrium


using ..Geometry


using Elliptic



include("MagneticField.jl")

export BFieldT
export compute_B!
export compute_B_isl!



include("qProfiles.jl")

export fu_dam_q
export qfm_benchmark_q



include("DensityProfiles.jl")

export uniform_dens


end
