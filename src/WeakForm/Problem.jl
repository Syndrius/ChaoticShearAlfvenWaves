"""
    init_problem(; fields::FieldsT, geometry::GeometryT, flr::FLRT=ideal_flr)

Initialises the problem struct, the key struct that dicates the problem to be solved.
"""
function init_problem(; fields::FieldsT, geometry::GeometryT, flr::FLRT=ideal_flr)
    
    return inst_problem(fields, geometry, flr) 
end

"""
    inst_problem(fields::IslandFieldsT, geo::GeometryT, flr::FLRT)

Instantiates the problem for island coordinates.
"""
function inst_problem(fields::IslandFieldsT, geo::GeometryT, flr::FLRT)
    inst_q(κ) = island_q(κ, fields.isls[1])
    inst_fields = IslandFieldsT(island_q, inst_q, fields.dens, fields.isls)
    inst_met(met, κ, ᾱ, τ, R0) = island_metric!(met, κ, ᾱ, τ, R0, fields.isls[1])
    inst_geo = GeometryT(:κ, geo.R0, inst_met, geo.a, geo.B0)
    return ProblemT(inst_fields, inst_geo, flr)
end


"""
    inst_problem(fields::RadialFieldsT, geo::GeometryT, flr::FLRT)

Instantiates the problem for geometric radius.
"""
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

    #this q-profile uses the island.
    if fields.q == island_q
        inst_q(r) = island_q(r, fields.isl[1])
        inst_fields = FieldsT()
        #metric is fine.
        return ProblemT(inst_fields, inst_geo, prob.flr)
    end

    return ProblemT(fields, inst_geo, flr)
end


"""
    inst_problem(fields::FluxFieldsT, geo::GeometryT, flr::FLRT)

Instantiates the problem for toroidal flux.
"""
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
    
    #this q-profile uses the island.
    if fields.q == island_q
        inst_q(ψ) = island_q(ψ, fields.isl[1])
        inst_fields = FieldsT()
        return ProblemT(inst_fields, inst_geo, prob.flr)
    end

    return ProblemT(fields, inst_geo, flr)
end
