
"""
Struct storing the geometrical parameters.

### Fields
- R0::Float64 - The major radius.
- a::Float64=1.0 - The minor radius. Not implemented yet, assumed to be 1.0 everywhere.
- B0::Float64=1.0 - The magnetic field strength at the axis. Not implemented yet, assumed to be 1.0 everywhere.
"""
@kwdef struct GeoParamsT
    R0 :: Float64
    a :: Float64 = 1.0 #not implemented yet, assume 1.
    B0 :: Float64 = 1.0 #not implemented yet, assume 1.
end

"""
    init_geo(; R0::Float64)

Function to initialise the GeoParamsT struct. Currently very barebones.
"""
function init_geo(; R0::Float64)

    return GeoParamsT(R0=R0)
end
