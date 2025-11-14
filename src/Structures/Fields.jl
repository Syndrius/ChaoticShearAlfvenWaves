
abstract type FieldsT end

#hopefully, we don't need 3 types of islands anymore!

struct FluxFieldsT <: FieldsT
    q_func :: Function #we may additionally want the option of q to be passed as a polynomail
    q :: FunctionWrapper{Tuple{Float64, Float64}, Tuple{Float64}}
    dens :: FunctionWrapper{Float64, Tuple{Float64}} #all dens functions follow the pattern.
    isls :: Array{FluxIslandT} 
end

struct IslandFieldsT <: FieldsT
    q_func :: Function #we may additionally want the option of q to be passed as a polynomail
    q :: FunctionWrapper{Tuple{Float64, Float64}, Tuple{Float64}}
    dens :: FunctionWrapper{Float64, Tuple{Float64}} #all dens functions follow the pattern.
    isls :: Array{CoordIslandT} 
end

struct RadialFieldsT <: FieldsT
    q_func :: Function #we may additionally want the option of q to be passed as a polynomail
    q :: FunctionWrapper{Tuple{Float64, Float64}, Tuple{Float64}}
    dens :: FunctionWrapper{Float64, Tuple{Float64}} #all dens functions follow the pattern.
    isls :: Array{RadialIslandT} 
end


#=
"""
Struct for storing the magnetic field and related variables at a given coordinate.

### Fields
- B::Array{Float64} - Vector storing the magnetic field.
- b::Array{Float64} - Vector storing the normalised magnetic field.
- dB::Array{Float64, 2} V- ector storing the derivative of the magnetic field, second index refers to derivative coordinate.
- db::Array{Float64, 2} - Vector storing the derivative of the normalised magnetic field, second index refers to derivative coordinate.
- mag_B::Array{Float64} - Magnitude of the magnetic field. Stored as an array so struct is immutable.
- dmag_B::Array{Float64} - Derivative of the magnitude of B, index refers to derivative coordinate.
"""
struct BFieldT
    B :: vec3
    b :: vec3
    dB :: mat33
    db :: mat33
    mag_B :: vec1
    dmag_B :: vec3
    function BFieldT()
        B = zeros(vec3)
        b = zeros(vec3)
        dB = zeros(mat33)
        db = zeros(mat33)
        mag_B = zeros(vec1)
        dmag_B = zeros(vec3)
        new(B, b, dB, db, mag_B, dmag_B)
    end
end
=#
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


#this will have to be moved within fields.
#for us to use defaults.
#=
function init_fields(; q::Function, dens::Function, isls::Array{<:IslandT})

    #can't really remember how the fkn function wrapper shite works.
    return FluxFieldsT(q, q, dens, dens, isls)
end
=#
