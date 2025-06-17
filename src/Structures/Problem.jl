#still a bit unsure if we need all these different problemos, 
#currenlt it is required for the array of islands to work
#but that seems like it is possible to work around
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
@kwdef struct RadProblemT <: ProblemT
    q_func :: Function #we may additionally want the option of q to be passed as a polynomail
    q :: FunctionWrapper{Tuple{Float64, Float64}, Tuple{Float64}}
    met_func :: Function = rad_toroidal_metric!
    met :: FunctionWrapper{Nothing, Tuple{MetT, Float64, Float64, Float64, Float64}}
    dens_func :: Function = uniform_dens #this should be done in the same way as q
    dens :: FunctionWrapper{Float64, Tuple{Float64}}
    isls :: Array{RadIslandT} = RadIslandT[] # very unfor that this doesn't work
    geo :: GeoParamsT
    flr :: FLRT = no_flr
end


@kwdef struct FluxProblemT <: ProblemT
    q_func :: Function #we may additionally want the option of q to be passed as a polynomail
    q :: FunctionWrapper{Tuple{Float64, Float64}, Tuple{Float64}}
    met_func :: Function = rad_toroidal_metric!
    met :: FunctionWrapper{Nothing, Tuple{MetT, Float64, Float64, Float64, Float64}}
    dens_func :: Function = uniform_dens #this should be done in the same way as q
    dens :: FunctionWrapper{Float64, Tuple{Float64}}
    isls :: Array{FluxIslandT} = FluxIslandT[] # very unfor that this doesn't work
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
    q_func :: Function #we may additionally want the option of q to be passed as a polynomail
    q :: FunctionWrapper{Tuple{Float64, Float64}, Tuple{Float64}}
    met_func :: Function = rad_toroidal_metric!
    met :: FunctionWrapper{Nothing, Tuple{MetT, Float64, Float64, Float64, Float64}}
    dens_func :: Function = uniform_dens #this should be done in the same way as q
    dens :: FunctionWrapper{Float64, Tuple{Float64}}
    isls :: Array{CoordIslandT} = CoordIslandT[] # very unfor that this doesn't work
    geo :: GeoParamsT
    flr :: FLRT = no_flr
end
#think we may not actually need this anymore.



#constant island storing the case without an island.
const no_rad_isl = RadIslandT(m0=1.0, n0=1.0, A=0.0)
const no_flux_isl = FluxIslandT(m0=1.0, n0=1.0, A=0.0)
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
function init_problem(; type::Symbol=:radial, q::Function, met::Symbol=:torus, dens::Function=uniform_dens, isl::IslandT=no_rad_isl, isls::Array{IslandT}=IslandT[], geo::GeoParamsT, flr::FLRT=no_flr)

    #this is fked beyond beleif.
    #this function is a fkn disaster again. need to have a flux boolean I think!
    #need to actually just change how this works, ideally so we can instantiate the problem later!
    #in particular the islands and the q-profiles.
    #puts island into array for consistent usage throughout code.
    if isempty(isls)
        isls = [isl]
    end

    #think ideally we would be open to having met as a funtion not just a symbol
    #in case anyone wants to create there own metric function
    if type == :radial

        new_isls = RadIslandT[]
        for isl in isls
            if isl isa FluxIslandT
                isl = convert_isl(isl)
            end
            push!(new_isls, isl)
        end

        if met==:torus
            met_func = rad_toroidal_metric!
        elseif met == :cylinder
            met_func = rad_cylindrical_metric!
        elseif met == :slab
            met_func = slab_metric! #same for both cases?
        else
            display("Metric not recognised.")
        end
        un_inst_prob = RadProblemT(q_func=q, met_func = met_func, dens_func = dens, q=q_profile, met=metric!, dens=density_profile, isls=new_isls, geo=geo, flr=flr)

    elseif type == :flux
        
        new_isls = FluxIslandT[]
        for isl in isls
            if isl isa RadIslandT
                isl = convert_isl(isl)
            end
            push!(new_isls, isl)
        end


        if met==:torus
            met_func = flux_toroidal_metric!
        elseif met == :cylinder
            met_func = flux_cylindrical_metric!
        elseif met == :slab
            met_func = slab_metric! #same for both cases?
        else
            display("Metric not recognised.")
        end
        un_inst_prob = FluxProblemT(q_func=q, met_func = met_func, dens_func = dens, q=q_profile, met=metric!, dens=density_profile, isls=new_isls, geo=geo, flr=flr)

    elseif type == :island
        #this should probably assert the q-profile is correct
        #and the islands are correct etc
        met_func = island_metric!

        new_isls = CoordIslandT[]
        for isl in isls
            if isl isa FluxIslandT
                isl = convert_isl(isl)
            elseif isl isa RadIslandT
                #this probably won't actually do anything!
                isl = convert_isl(isl)
            end
            push!(new_isls, isl)
        end

        un_inst_prob = IslProblemT(q_func=q, met_func = met_func, dens_func = dens, q=q_profile, met=metric!, dens=density_profile, isls=new_isls, geo=geo, flr=flr)

    else
        display("Only radial, flux or island problem are currently implemented.")

    end


    return inst_problem(un_inst_prob)
end


#=
    #ok so ideally we would want to be ablet to initialise the q-profile and the island metric
    #this will require some fixes,

    #think we should just have a cylindrical met, and a toroidal met
    #but they are multiple dispatch based on the island being rad or flux
    #This should be easy enough to implement, as we can define the islands as flux or not, 
    #and the no_isl can be a flux or rad.
    #then we can do the same thing with B.
    #Not sure if we actually want the option to compute B in other weird ways
    #we can probably also have an IslandCoords islant to keep up the multiple dispatch
    #although that does still need the island metric so who knows.
    #unless we want to make a anon function for the island?
    #which will require getting anon-functions to work properly.
    #perhaps we can just make the problem store 2 q profiles etc?
    #not an ideal situation
    #but will allow the original q profile to be written then reconstructed as needed.

    #ideally we will do a better job with qfm in the future. Ideally, there would only be a single construct
    #for each grid type.
    #I think in the long term future, different compute B functions would be specified in the metric
    #it is somewhat restrictive atm, assumes certain properties of the coordinates which may not always be true
    if met == :torus
        met_func = rad_toroidal_metric!
    elseif met == :cylinder
        #are we doing radius? Who knows...
        met_func = rad_cylindrical_metric!
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
        #this is beyond useless now, JLD2 cannot save anon functions.
        isl_q_prof(κ::Float64) = island_coords_q(κ, isl)
        #return IslProblemT(q=isl_q_prof, met=island_metric!, dens=dens, isls=[isl], flr=flr, geo=geo)
        #we have now just hardcoded the island q-profile for specific islands,
        #as jld2 cannot write the anon-function
        #this is an awful solution 
        #TODO
        return IslProblemT(q=q, met=island_metric!, dens=dens, isls=[isl], flr=flr, geo=geo, B=isl_compute_B!)
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
        return TorProblemT(q=q_prof, met=met_func, dens=dens, isls=[isl], flr=flr, geo=geo, B=B)
    else
        new_isls = []
        if isls[1] != no_isl || length(isls) > 1
            #arguably don't wont to do this for everycase
            #see qfm benchmark for example when this is annoying af.
            #may need a try catch or something
            try
                for i in 1:length(isls)
                    #isls[i] = inst_island(isls[i], q)
                    push!(new_isls, inst_island(isls[i], q))
                end
            catch e
                new_isls = isls
                display(new_isls)
                display("Unable to instantiate islands.")
            end

        end
        return TorProblemT(q=q, met=met_func, dens=dens, isls=new_isls, flr=flr, geo=geo, B = B)
    end 

end
=#

function inst_problem(prob::IslProblemT)

    if length(prob.isls) > 1
        display("Island q-profile only supports a single island.")
        return
    end

    #not sure if this actually happens for all of these.
    isl = inst_island(prob.isls[1])
    if isnan(isl.w) || isnan(isl.A)
        display("Please specify m0, n0, r0 and either w or A.")
        return
    end
    inst_island_q_prof(r::Float64) = island_coords_q(r, isl)
    #probably should assert this matches the proper q-profile!
    
    #not sure this knows about the above instantiated form of the q-profile.
    inst_island_metric(met::MetT, κ::Float64, ᾱ::Float64, τ::Float64, R0::Float64) = island_metric!(met, κ, ᾱ, τ, R0, isl)
    return IslProblemT(q_func=prob.q_func, met_func = prob.met_func, dens_func = prob.dens_func, q=inst_island_q_prof, met=inst_island_metric, dens=prob.dens_func, isls=[isl], geo=prob.geo, flr=prob.flr)

end

#actually creates a usable problemo
#mainly used for creating anon-function etc
#but also for saving/loading problem from file.
#this needs to be split for each problem type
function inst_problem(prob::RadProblemT)
    
    #note that this does not handle flux coordinates at all.
    #TODO
    if prob.q_func == island_equiv_q

        if length(prob.isls) > 1
            display("Island q-profile only supports a single island.")
            return
        end

        #think this will be the only case where we insantiate the island tbh!
        #we only really use the istantiated form when doing mapping
        #which we probably will only ever do with this q-profile.
        isl = inst_island(prob.isls[1])
        if isnan(isl.w) || isnan(isl.A)
            display("Please specify m0, n0, r0 and either w or A.")
            return
        end
        inst_equiv_q_prof(r::Float64) = island_equiv_q(r, isl)

        #the array of isl won't work, rip.
        return RadProblemT(q_func=prob.q_func, met_func = prob.met_func, dens_func = prob.dens_func, q=inst_equiv_q_prof, met=prob.met_func, dens=prob.dens_func, isls=[isl], geo=prob.geo, flr=prob.flr)

        #need some assertation about flux vs radial islands!
        #and or any other islands we use.
        #i.e. the islandcoords case.
        #note that both of these should be true or neither.
    elseif prob.q_func == island_coords_q || prob.met_func == island_metric!
        if length(prob.isls) > 1
            display("Island q-profile only supports a single island.")
            return
        end

        #not sure if this actually happens for all of these.
        isl = inst_island(prob.isls[1])
        if isnan(isl.w) || isnan(isl.A)
            display("Please specify m0, n0, r0 and either w or A.")
            return
        end
        inst_island_q_prof(r::Float64) = island_coords_q(r, isl)
        #probably should assert this matches the proper q-profile!
        
        #not sure this knows about the above instantiated form of the q-profile.
        inst_island_metric(met::MetT, κ::Float64, ᾱ::Float64, τ::Float64, R0::Float64) = island_metric!(met, κ, ᾱ, τ, R0, isl)
        return RadProblemT(q_func=prob.q_func, met_func = prob.met_func, dens_func = prob.dens_func, q=inst_island_q_prof, met=inst_island_metric, dens=prob.dens_func, isls=[isl], geo=prob.geo, flr=prob.flr)

    else
        return RadProblemT(q_func=prob.q_func, met_func = prob.met_func, dens_func = prob.dens_func, q=prob.q_func, met=prob.met_func, dens=prob.dens_func, isls=prob.isls, geo=prob.geo, flr=prob.flr)
    end

    #not doing much with the islands atm!
    #return ProblemT(q_func=prob.q_func, met_func = prob.met_func, dens_func = prob.dens_func, q=inst_q_prof, met=inst_metric, dens=inst_dens, isls=prob.isls, geo=prob.geo, flr=prob.flr)

end

#this is a disaster.
function inst_problem(prob::FluxProblemT)
    
    #note that this does not handle flux coordinates at all.
    #TODO
    #don't think this exists in the flux case
    if prob.q_func == island_equiv_q

        if length(isls) > 1
            display("Island q-profile only supports a single island.")
            return
        end

        #think this will be the only case where we insantiate the island tbh!
        #we only really use the istantiated form when doing mapping
        #which we probably will only ever do with this q-profile.
        isl = inst_island(isls[1])
        if isnan(isl.w) || isnan(isl.A)
            display("Please specify m0, n0, r0 and either w or A.")
            return
        end
        inst_equiv_q_prof(r::Float64) = island_equiv_q(r, isl)

        #the array of isl won't work, rip.
        return FluxProblemT(q_func=prob.q_func, met_func = prob.met_func, dens_func = prob.dens_func, q=ins_equiv_q_prof, met=prob.met_func, dens=prob.dens_func, isls=[isl], geo=prob.geo, flr=prob.flr)

        #need some assertation about flux vs radial islands!
        #and or any other islands we use.
        #i.e. the islandcoords case.
        #note that both of these should be true or neither.
        #pretty sure this doesn't work in flux either!
    elseif prob.q_func == island_coords_q || prob.met_func == island_metric!
        if length(prob.isls) > 1
            display("Island q-profile only supports a single island.")
            return
        end

        #not sure if this actually happens for all of these.
        isl = inst_island(prob.isls[1])
        if isnan(isl.w) || isnan(isl.A)
            display("Please specify m0, n0, r0 and either w or A.")
            return
        end
        inst_island_q_prof(r::Float64) = island_coords_q(r, isl)
        #probably should assert this matches the proper q-profile!
        
        #not sure this knows about the above instantiated form of the q-profile.
        inst_island_metric(met::MetT, κ::Float64, ᾱ::Float64, τ::Float64, R0::Float64) = island_metric!(met, κ, ᾱ, τ, R0, isl)
        return FluxProblemT(q_func=prob.q_func, met_func = prob.met_func, dens_func = prob.dens_func, q=inst_island_q_prof, met=inst_island_metric, dens=prob.dens_func, isls=[isl], geo=prob.geo, flr=prob.flr)

    else
        return FluxProblemT(q_func=prob.q_func, met_func = prob.met_func, dens_func = prob.dens_func, q=prob.q_func, met=prob.met_func, dens=prob.dens_func, isls=prob.isls, geo=prob.geo, flr=prob.flr)
    end

    #not doing much with the islands atm!
    #return ProblemT(q_func=prob.q_func, met_func = prob.met_func, dens_func = prob.dens_func, q=inst_q_prof, met=inst_metric, dens=inst_dens, isls=prob.isls, geo=prob.geo, flr=prob.flr)

end
