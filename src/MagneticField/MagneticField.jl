

module MagneticField

using MID.Geometry


include("ComputeField.jl")

export compute_B!
export BFieldT


include("qProfiles.jl")

export Axel_q


#maybe combine into just profiles?
include("DensityProfiles.jl") #weird spot!

export uniform_dens

end