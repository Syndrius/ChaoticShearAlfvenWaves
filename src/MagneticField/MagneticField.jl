"""

Contains the Magnetic field and functions to calculate it, includes
 - BFieldT struct and function to fill it.
 - Variety of q and density profiles that are passed in as arguments.

"""
module MagneticField

using MID.Geometry


include("ComputeField.jl")

export BFieldT
export init_empty_B
export compute_B!
export compute_island_B!


include("qProfiles.jl")

export Axel_q
export island_damping_q
export bowden_singular_q
export comparison_bowden_q
export fu_dam_q
export default_island_q
export island_3_2_q
export symmetric_q
export flr_q
export test_q
export Axel_island_q
export island_mode_q
export island_mode_21
export gae_isl_q


include("DensityProfiles.jl") #weird spot!

export uniform_dens
export axel_dens
export bowden_singular_dens
export comparison_bowden_dens
export gae_isl_dens

end