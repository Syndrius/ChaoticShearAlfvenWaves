


"""
Struct for storing the metric which describes the geometry and related derivatives at each coordinate.

### Fields
- gl::Array{Float64, 2} - 3x3 matrix storing the lowered metric g_{ij}
- gu::Array{Float64, 2} - 3x3 matrix storing the raised metric g^{ij}
- dgl::Array{Float64, 3} - Derivative of gl, 3rd index labels coordinate that derivative is taken with respect to.
- dgl::Array{Float64, 3} - Derivative of gu, 3rd index labels coordinate that derivative is taken with respect to.
- J::Array{Float64} - Jacobian of the metric. Stored as an array so struct is immutable.
- dJ::Array{Float64, 1} - Derivative of J, index labels coordinate that derivative is taken with respect to.
"""
#=
struct MetT
    gl :: mat33
    gu :: mat33
    dgl :: mat333
    dgu :: mat333
    J :: vec1
    dJ :: vec3
    function MetT()
        #kind of dumb
        #new(zeros(MMatrix{3, 3}), zeros(MMatrix{3, 3}), zeros(MArray{3, 3, 3}), zeros(MArray{3, 3, 3}), zeros(MVector{1}), zeros(MVector{3}))
        gl = zeros(mat33)
        gu = zeros(mat33)
        dgl = zeros(mat333)
        dgu = zeros(mat333)
        J = zeros(vec1)
        dJ = zeros(vec3)
        new(gl, gu, dgl, dgu, J, dJ)
    end
end
=#

struct MetT
    gl :: Array{Float64, 2}
    gu :: Array{Float64, 2} 
    dgl :: Array{Float64, 3} 
    dgu :: Array{Float64, 3} 
    J :: Array{Float64, 1}
    dJ :: Array{Float64, 1} 
    function MetT()
        new(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3), zeros(1), zeros(3))
    end
end


"""
Struct storing the geometrical parameters.

### Fields
- R0::Float64 - The major radius.
- a::Float64=1.0 - The minor radius. Not implemented yet, assumed to be 1.0 everywhere.
- B0::Float64=1.0 - The magnetic field strength at the axis. Not implemented yet, assumed to be 1.0 everywhere.
"""
struct GeometryT
    type :: Symbol
    R0 :: Float64
    met :: FunctionWrapper{Nothing, Tuple{MetT, Float64, Float64, Float64, Float64}}
    a :: Float64 #not implemented yet, assume 1.
    B0 :: Float64 #not implemented yet, assume 1.
end


#=
"""

If we want this to have default, this will need to be moved into geometry.
"""
function init_geometry(R0::Float64, met::Function)

    return GeometryT(R0, met, met, 1.0, 1.0)
end
=#
