"""
This struct stores Finite Larmor Radius effects to add.

### Fields
- δ::Float64=0.0 - Artifical damping.  
- ρ_i::Float64=0.0 - Ion gyro radius.
- δ_e::Float64=0.0 - Electron damping.
"""
@kwdef struct FLRT
    δ :: Float64 = 0.0 
    ρ_i :: Float64 = 0.0
    δ_e :: Float64 = 0.0
end


"""
    init_flr(;δ=0.0::Float64, ρ_i=0.0::Float64, δ_e=0.0::Float64)

Initialise struct for storing any additional finite Larmor radius effects.
Includes artificial resistivity, δ, ion radius, ρ_i, electron damping, δ_e.
"""
function init_flr(;δ=0.0::Float64, ρ_i=0.0::Float64, δ_e=0.0::Float64)
    return FLRT(δ, ρ_i, δ_e)
end

