

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
This struct stores the key information that defines the problem to be solved.
One of the main inputs for matrix construction functions.

### Fields
- q::Function - The q-profile, function of r that return (q, dq)
- met::Function=toroidal_metric! - Computes the metric at each coordinates, defaults to toroidal.
- dens::Function=uniform_dens - Computes the density as a function of r, defaults to uniform density.
- isl::IslandT=no_isl Struct - storing the island parameters, defaults to no island.
- geo::GeoParamsT - Struct storing the geometrical parameters.
- flr::FLRT - Struct storing finite larmor effects.
"""
@kwdef struct TorProblemT <: ProblemT
    q :: Function 
    compute_met :: Function = toroidal_metric!
    dens :: Function = uniform_dens
    isl :: IslandT = no_isl
    isl2 :: IslandT = no_isl
    geo :: GeoParamsT
    flr :: FLRT = no_flr
end


#problem with (κ, ᾱ, φ)
@kwdef struct IslProblemT <: ProblemT
    dens :: Function = uniform_dens
    geo :: GeoParamsT
    flr :: FLRT = no_flr
    isl :: IslandT #need at least m,n, r0 and A/w.
end



#constant island storing the case without an island.
const no_isl = IslandT(m0=1.0, n0=1.0, A=0.0)
#constant flr for cases without any flr corrections.
const no_flr = FLRT(δ=0.0, ρ_i=0.0, δ_e=0.0)

"""
Constructor for struct which stores the key information that defines the problem to be solved.
Main input for matrix construction functions.

# Args
- q::Function The q-profile, function of r that return (q, dq)
- met::Function=toroidal_metric! Computes the metric at each coordinates, defaults to toroidal.
- dens::Function=uniform_dens Computes the density as a function of r, defaults to uniform density.
- isl::IslandT=no_isl Struct storing the island parameters, defaults to no island.
- geo::GeoParamsT Struct storing geometrical parameters.
- flr::FLRT=no_flr - Struct storing finite larmor effects, defaults to no corrections.. 
"""
function init_problem(; q::Function, met::Function=toroidal_metric!, dens::Function=uniform_dens, isl::IslandT=no_isl, isl2::IslandT=no_isl,geo::GeoParamsT, flr::FLRT=no_flr)

    #may want to indroduce symbold into here
    #for the metric :toroidal etc may be clearer.
    if length(methods(q)[1].sig.parameters) == 3
        #this q -profile uses island parameters as input, this creates an appropriate q-profile
        #we should probably be asserting some stuff about the ol island first
        q_prof(r) = q(r, isl)
        #not sure if this actually works
        return ProblemT(q=q_prof, compute_met=met, dens=dens, isl=isl, flr=flr, geo=geo)
    else
        if isl != no_isl
            isl = inst_island(isl, q)
        end
        if isl2 != no_isl
            isl2 = inst_island(isl2, q)
        end
        return TorProblemT(q=q, compute_met=met, dens=dens, isl=isl, isl2=isl2, flr=flr, geo=geo)
    end 

    
end


function init_isl_problem(; dens::Function=uniform_dens, geo::GeoParamsT, flr::FLRT=no_flr, isl::IslandT)

    isl = inst_island(isl)

    return IslProblemT(dens=dens, geo=geo, flr=flr, isl=isl)

    
end