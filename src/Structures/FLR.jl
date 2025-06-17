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

#TODO
function init_flr()
end

