"""

Module for the equilibrium state of the plasma. This includes the density and the magnetic field, computed from an input q-profile and input magnetic islands.
"""
module Equilibrium


using ..Geometry


using Elliptic
using FunctionWrappers
import FunctionWrappers: FunctionWrapper



include("MagneticField.jl")

export BFieldT
export compute_B!
#export compute_B_isl!



include("qProfiles.jl") #once again, this has become a disaster. Ideally, we should just have a few q-profile written, and then perhaps we read a file for extras or something

export q_profile
export fu_dam_q
export low_shear_q
export low_shear_qfm_q
export qfm_q
export qfm_benchmark_q
export island_q



include("DensityProfiles.jl")

export density_profile
export uniform_dens


end
