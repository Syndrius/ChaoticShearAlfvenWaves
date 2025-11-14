
#should this struct contain a field for Hermitian?
#and should we bother with the SAWProblem anymore?
#probably not.
function init_problem(; fields::FieldsT, geometry::GeometryT, flr::FLRT=ideal_flr)
    #something like
    #do we just not deal with this.
    #let the user fail.
    #storing a var is kind of fkn stupid.
    #just not going to bother about this I think!
    #we could revert back to an island_W_and_I!, but this is not ideal for the q-profile...
    #although
    #this is such a fkn awkward problem, really just for the single cooked af q-profile
    return inst_problem(fields, geometry, flr) #how is post-processing going to have access to this!
end

function inst_problem(fields::IslandFieldsT, geo::GeometryT, flr::FLRT)
    inst_q(κ) = island_q(κ, fields.isls[1])
    inst_fields = IslandFieldsT(island_q, inst_q, fields.dens, field.isls)
    inst_met(met, κ, ᾱ, τ) = island_metric!(met, κ, ᾱ, τ, fields.isls[1])
    inst_geo = GeometryT(:κ, geo.R0, inst_met, geo.a, geo.B0)
    return ProblemT(inst_fields, inst_geo, prob.flr)
end


#NOTE!
#don't think we actually have two slab metrics
#might just remove tbh!
function inst_problem(fields::RadialFieldsT, geo::GeometryT, flr::FLRT)

    #this is done here to cover the case where metric depends on the fields
    #eg for island case.
    if geo.type == :cyl
        inst_geo = GeometryT(:cyl, geo.R0, radial_cylindrical_metric!, geo.a, geo.B0)
    elseif geo.type == :cart
        inst_geo = GeometryT(:cart, geo.R0, radial_cartesian_metric!, geo.a, geo.B0)
    else
        inst_geo = GeometryT(:tor, geo.R0, radial_toroidal_metric!, geo.a, geo.B0)
    end
    

    if fields.q == island_q
        #may need a check on the island here tbh!
        #and a computation if any are needed
        #think this assumes that fields checks that the island is in the same form
        #and the magnetic field. -> but maybe not.
        #want to assert only a single island is used.
        inst_q(r) = island_q(r, fields.isl[1])
        inst_fields = FieldsT()
        #metric is fine.
        return ProblemT(inst_fields, inst_geo, prob.flr)
    end

    #this is pretty fkn awkward.
    return ProblemT(fields, inst_geo, flr)
end


function inst_problem(fields::FluxFieldsT, geo::GeometryT, flr::FLRT)
    #this is done here to cover the case where metric depends on the fields
    #eg for island case.
    if geo.type == :cyl
        inst_geo = GeometryT(:cyl, geo.R0, flux_cylindrical_metric!, geo.a, geo.B0)
    elseif geo.type == :cart
        inst_geo = GeometryT(:cart, geo.R0, flux_cartesian_metric!, geo.a, geo.B0)
    else
        inst_geo = GeometryT(:tor, geo.R0, flux_toroidal_metric!, geo.a, geo.B0)
    end
    

    if fields.q == island_q
        #may need a check on the island here tbh!
        #and a computation if any are needed
        #think this assumes that fields checks that the island is in the same form
        #and the magnetic field. -> but maybe not.
        #want to assert only a single island is used.
        inst_q(ψ) = island_q(ψ, fields.isl[1])
        inst_fields = FieldsT()
        #metric is fine.
        return ProblemT(inst_fields, inst_geo, prob.flr)
    end

    #this is pretty fkn awkward.
    return ProblemT(fields, inst_geo, flr)
end
