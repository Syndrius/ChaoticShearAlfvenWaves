"""

Contains the Magnetic field and functions to calculate it, includes
 - BFieldT struct and function to fill it.
 - Variety of q and density profiles that are passed in as arguments.

"""
module Equilibrium

using ..Geometry

using Elliptic



include("MagneticField.jl")

export BFieldT
export compute_B!
export compute_B_isl!



include("qProfiles.jl")

export Axel_q
#export island_damping_q
#export bowden_singular_q
#export comparison_bowden_q
export fu_dam_q
export qfm_benchmark_q
#export default_island_q
#export island_3_2_q
#export symmetric_q
#export flr_q
#export test_q
#export Axel_island_q
#export island_mode_q
#export island_mode_21
#export inside_island_q
#export gae_isl_q
#export tae_isl_damping_q
#export chaos_q
#export qfm_q



include("DensityProfiles.jl")

export uniform_dens
#export axel_dens
#export bowden_singular_dens
#export comparison_bowden_dens
#export gae_isl_dens

end
