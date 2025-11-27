#Basic tests with default inputs.
#should define the problem here!

geo = init_geometry()
fields = init_fields()


prob = init_problem(fields=fields, geometry=geo)

display("Running Flux Continuum Test")
include("FluxContinuum.jl")

display("Running Flux FSS Test")
include("FluxFSS.jl")

display("Running Flux FFS Test")
include("FluxFFS.jl")

display("Running Flux FFF Test")
include("FluxFFF.jl")

##################################

geo = init_geometry()
fields = init_fields(:r)

prob = init_problem(geometry=geo, fields=fields)

display("Running Radial Continuum Test")
include("RadialContinuum.jl")

display("Running Radial FSS Test")
include("RadialFSS.jl")

display("Running Radial FFS Test")
include("RadialFFS.jl")

display("Running Radial FFF Test")
include("RadialFFF.jl")
