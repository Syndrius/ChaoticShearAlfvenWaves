

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
Struct storing the geometry.

### Fields
- type::Symbol - denotes the type of geometry, :tor, :cyl or :isl.
- R0::Float64 - The major radius.
- met :: FunctionWrapper{Nothing, Tuple{MetT, Float64, Float64, Float64, Float64}} - function for computing the metric.
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

