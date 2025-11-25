#abstract type that stores the fields used 
#note that dipatch is performed on island type, this is used for future proofing.
abstract type FieldsT end

"""
Three differenent cases for different radial coordinate used.
Note that the q-function is defined twice for creating anonymous functions that can be written to file.

### Fields
- q_func::Function - Generic function, stored for instantiation.
- q::FunctionWrapper{Tuple{Float64, Float64}, Tuple{Float64}} - Actual function used in computation.
- dens::FunctionWrapper{Float64, Tuple{Float64}} - Function to compute the density.
- isls::Array{FluxIslandT} - Islands used, the type is different for multiple dispatch.


"""

struct FluxFieldsT <: FieldsT
    q_func :: Function 
    q :: FunctionWrapper{Tuple{Float64, Float64}, Tuple{Float64}}
    dens :: FunctionWrapper{Float64, Tuple{Float64}} 
    isls :: Array{FluxIslandT} 
end

struct IslandFieldsT <: FieldsT
    q_func :: Function 
    q :: FunctionWrapper{Tuple{Float64, Float64}, Tuple{Float64}}
    dens :: FunctionWrapper{Float64, Tuple{Float64}} 
    isls :: Array{CoordIslandT} 
end

struct RadialFieldsT <: FieldsT
    q_func :: Function 
    q :: FunctionWrapper{Tuple{Float64, Float64}, Tuple{Float64}}
    dens :: FunctionWrapper{Float64, Tuple{Float64}} 
    isls :: Array{RadialIslandT} 
end


"""
Struct for storing the magnetic field and related variables at a given coordinate.

### Fields
- B::Array{Float64} - Vector storing the magnetic field.
- b::Array{Float64} - Vector storing the normalised magnetic field.
- dB::Array{Float64, 2} - Vector storing the derivative of the magnetic field, second index refers to derivative coordinate.
- db::Array{Float64, 2} - Vector storing the derivative of the normalised magnetic field, second index refers to derivative coordinate.
- mag_B::Array{Float64} - Magnitude of the magnetic field. Stored as an array so struct is immutable.
- dmag_B::Array{Float64} - Derivative of the magnitude of B, index refers to derivative coordinate.
"""
struct BFieldT
    B :: Array{Float64, 1} 
    b :: Array{Float64, 1} 
    dB :: Array{Float64, 2} 
    db :: Array{Float64, 2} 
    mag_B :: Array{Float64, 1}
    dmag_B :: Array{Float64, 1} 
    function BFieldT()
        new(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), zeros(1), zeros(3))
    end
end

