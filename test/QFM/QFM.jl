
geo = init_geometry(:tor, R0=4.0)

#probably cant just use w huh.
#this would require the root solve shite.
isl = init_island(m0=3, n0=1, A=0.1)

fields = init_fields(:ψ, q=cantori_q, isl=isl)

prob = init_problem(geometry=geo, fields=fields)

dir = abspath(joinpath(pathof(ChaoticShearAlfvenWaves), "../../test/data/"))
surfs = surfaces_from_file(joinpath(dir, "benchmark_surfaces.jld2"))


display("Running Continuum Test")
include("Continuum.jl")

display("Running FSS Test")
include("FSS.jl")

display("Running FFS Test")
include("FFS.jl")

display("Running FFF Test")
include("FFF.jl")

#perhaps a mapping test?
