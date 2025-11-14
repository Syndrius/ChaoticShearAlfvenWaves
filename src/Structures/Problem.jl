#this is a fookin disaster.

#probably don't need this!

struct ProblemT 
    fields :: FieldsT
    geo :: GeometryT
    flr :: FLRT
end
#still a bit unsure if we need all these different problemos, 
#currenlt it is required for the array of islands to work

#constant island storing the case without an island.
const no_rad_isl = RadialIslandT(m0=1.0, n0=1.0, A=0.0, r0=1.0, qp=1.0, q0=1.0, w=0.0)
const no_flux_isl = FluxIslandT(m0=1.0, n0=1.0, A=0.0) #unused I think!
#constant flr for cases without any flr corrections.
const ideal_flr = FLRT(δ=0.0, ρ_i=0.0, δ_e=0.0)

#function init_problem(; fields::FieldsT, geo::GeometryT, flr::FLRT=ideal_flr)
#    return SAWProblemT(fields, geo, flr)
#end
#the whole inst problem thingo will be annoying af.

