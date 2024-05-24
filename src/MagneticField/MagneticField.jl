

module MagneticField

using MID.Geometry


include("ComputeField.jl")

export compute_B!
export BFieldT


include("qProfiles.jl")

export Axel_q
export island_damping_q
export singular_bowden_q
export comparison_bowden_q


#maybe combine into just profiles?
include("DensityProfiles.jl") #weird spot!

export uniform_dens
export axel_dens
export bowden_singular_dens
export comparison_bowden_dens

end