"""
This struct stores Finite Larmor Radius effects to add.

### Fields
- δ::Float64 - Artifical damping.  
- ρ_i::Float64 - Ion gyro radius.
- δ_e::Float64 - Electron damping.
"""
struct FLRT
    δ :: Float64
    ρ_i :: Float64
    δ_e :: Float64
end

#flr for ideal MHD.
const ideal_flr = FLRT(0.0, 0.0, 0.0)

"""
    init_flr(;δ=0.0::Float64, ρ_i=0.0::Float64, δ_e=0.0::Float64)

Initialise struct for storing any additional finite Larmor radius effects.
Includes artificial resistivity, δ, ion radius, ρ_i, electron damping, δ_e.
"""
function init_flr(;δ=0.0::Float64, ρ_i=0.0::Float64, δ_e=0.0::Float64)
    return FLRT(δ, ρ_i, δ_e)
end

