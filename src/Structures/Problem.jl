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
- isl::IslandT=no_isl - Struct storing the island parameters, defaults to no island.
- geo::GeoParamsT - Struct storing the geometrical parameters.
- flr::FLRT - Struct storing finite larmor effects.
"""
@kwdef struct TorProblemT <: ProblemT
    q :: Function 
    met :: Function = toroidal_metric!
    dens :: Function = uniform_dens
    isls :: Array{IslandT} = IslandT[]
    geo :: GeoParamsT
    flr :: FLRT = no_flr
end


"""
Struct for problems using island coordinates, defined for multiple dispatch. Less information is needed here as the metric and q-profile are predetermined.

### Fields
- dens::Function=uniform_dens - Computes the density as a function of r, defaults to uniform density.
- geo::GeoParamsT - Struct storing the geometrical parameters.
- flr::FLRT - Struct storing finite larmor effects.
- isl::IslandT - Struct storing the island parameters, defaults to no island.
"""
@kwdef struct IslProblemT <: ProblemT
    q :: Function # This will automatically be the island q
    met :: Function = island_metric!
    dens :: Function = uniform_dens
    geo :: GeoParamsT
    flr :: FLRT = no_flr
    isls :: Array{IslandT} #need at least m,n, r0 and A/w., should always be a single island.
end



#constant island storing the case without an island.
const no_isl = IslandT(m0=1.0, n0=1.0, A=0.0)
#constant flr for cases without any flr corrections.
const no_flr = FLRT(δ=0.0, ρ_i=0.0, δ_e=0.0)


"""
Constructor for struct which stores the key information that defines the problem to be solved.
Main input for matrix construction functions.

# Args
- q::Function - The q-profile, function of r that return (q, dq)
- met::Function=toroidal_metric! - Computes the metric at each coordinates, defaults to toroidal.
- dens::Function=uniform_dens - Computes the density as a function of r, defaults to uniform density.
- isl::IslandT=no_isl - Struct storing the island parameters, defaults to no island.
- geo::GeoParamsT - Struct storing geometrical parameters.
- flr::FLRT=no_flr - Struct storing finite larmor effects, defaults to no corrections.. 
"""
function init_problem(; q::Function, met::Symbol=:torus, dens::Function=uniform_dens, isl::IslandT=no_isl, isls::Array{IslandT}=IslandT[], geo::GeoParamsT, flr::FLRT=no_flr)

    #puts island into array for consistent usage throughout code.
    if isempty(isls)
        isls = [isl]
    end

    #ideally we will do a better job with qfm in the future. Ideally, there would only be a single construct
    #for each grid type.
    #I think in the long term future, different compute B functions would be specified in the metric
    #it is somewhat restrictive atm, assumes certain properties of the coordinates which may not always be true
    if met == :torus
        met_func = toroidal_metric!
    elseif met == :cylinder
        met_func = cylindrical_metric!
        #bit annoying that we have to specify the q-profile, we may want to put a default one of fu-dam.
    elseif met == :island || q == island_coords_q
        if length(isls) > 1
            display("Island metric only supports a single island.")
            return
        end
        isl = inst_island(isls[1])
        if isnan(isl.w) || isnan(isl.A)
            display("Please specify m0, n0, r0 and either w or A.")
            return
        end
        isl_q_prof(κ::Float64) = island_coords_q(κ, isl)
        return IslProblemT(q=isl_q_prof, met=island_metric!, dens=dens, isls=[isl], flr=flr, geo=geo)
    elseif met == :slab #this may not be possible tbh!
        met_func = slab_metric!
    else
        display("Metric not available")
        return
    end

    #may want to indroduce symbold into here
    #for the metric :toroidal etc may be clearer.
    #TODO, no longer works now that islands are done via array.
    #no idea what the fek is going on here tbh!
    if length(methods(q)[1].sig.parameters) == 3 || q == island_equiv_q #not sure the first check actually works.
        if length(isls) > 1
            display("Constructing island q-profile only supports a single island.")
            return
        end
        #this is the case where we are using a q-profile that depends on the island
        #here it is converted to the standard input.
        #this is wrong tho!
        isl = inst_island(isls[1])
        if isnan(isl.w) || isnan(isl.A)
            display("Please specify m0, n0, r0 and either w or A.")
            return
        end
        display("are we here")
        q_prof(r::Float64) = island_equiv_q(r, isl)
        #not sure if this actually works
        return TorProblemT(q=q_prof, met=met_func, dens=dens, isls=[isl], flr=flr, geo=geo)
    else
        if isls[1] != no_isl || length(isls) > 1
            #arguably don't wont to do this for everycase
            #see qfm benchmark for example when this is annoying af.
            #may need a try catch or something
            new_isls = []
            try
                for i in 1:length(isls)
                    #isls[i] = inst_island(isls[i], q)
                    push!(new_isls, inst_island[i], q)
                end
            catch e
                new_isls = isls
                display(new_isls)
                display("Unable to instantiate islands.")
            end

        end
        return TorProblemT(q=q, met=met_func, dens=dens, isls=new_isls, flr=flr, geo=geo)
    end 

end


"""
    init_isl_problem(; dens::Function=uniform_dens, geo::GeoParamsT, flr::FLRT=no_flr, isl::IslandT)

Function to initialise an island problem.
"""
function init_isl_problem(; dens::Function=uniform_dens, geo::GeoParamsT, flr::FLRT=no_flr, isl::IslandT)

    isl = inst_island(isl)

    return IslProblemT(dens=dens, geo=geo, flr=flr, isl=isl)

    
end


"""
    init_geo(; R0::Float64)

Function to initialise the GeoParamsT struct. Currently very barebones.
"""
function init_geo(; R0::Float64)

    return GeoParamsT(R0=R0)
end
