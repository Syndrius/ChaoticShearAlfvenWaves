

"""
Struct storing the geometrical parameters.
# Fields
R0::Float64  - The major radius
a::Float64=1 - The minor radius
B0::Float64=1 - The magnetic field strength at the axis.
"""
@kwdef struct GeoParamsT
    R0 :: Float64
    a :: Float64 = 1.0 #not implemented yet, assume 1.
    B0 :: Float64 = 1.0 #not implemented yet, assume 1.
end


"""
This struct stores the key information that defines the problem to be solved.
Main input for matrix construction functions.

# Fields
q::Function - The q-profile, function of r that return (q, dq)
met::Function=toroidal_metric! - Computes the metric at each coordinates, defaults to toroidal.
dens::Function=uniform_dens - Computes the density as a function of r, defaults to uniform density.
isl::IslandT=no_isl Struct - storing the island parameters, defaults to no island.
geo::GeoParamsT - Struct storing the geometrical parameters.
δ::Float64=0.0 - Artifical damping.    
"""
@kwdef struct ProblemT
    q :: Function 
    compute_met :: Function = toroidal_metric!
    dens :: Function = uniform_dens
    isl :: IslandT = no_isl
    geo :: GeoParamsT
    δ :: Float64 = 0.0 # may be better to have this as a separate arg for construct if it changing a lot.
end





#constant island storing the case without an island.
const no_isl = IslandT(m0=1.0, n0=1.0, A=0.0)

"""
Constructor for struct which stores the key information that defines the problem to be solved.
Main input for matrix construction functions.

# Args
- q::Function The q-profile, function of r that return (q, dq)
- met::Function=toroidal_metric! Computes the metric at each coordinates, defaults to toroidal.
- dens::Function=uniform_dens Computes the density as a function of r, defaults to uniform density.
- isl::IslandT=no_isl Struct storing the island parameters, defaults to no island.
- R0::Float64 Major radius.
- δ::Float64=0.0 Artifical damping.    
"""
function init_problem(; q::Function, met::Function=toroidal_metric!, dens::Function=uniform_dens, isl::IslandT=no_isl, geo::GeoParamsT, δ::Float64=0.0)


    return ProblemT(q=q, compute_met=met, dens=dens, isl=isl, δ=δ, geo=geo)
end