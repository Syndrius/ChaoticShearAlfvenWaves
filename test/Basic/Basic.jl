#Basic tests with default inputs.

display("Running Flux Continuum Test")
include("FluxContinuum.jl")

display("Running Flux FSS Test")
include("FluxFSS.jl")

display("Running Flux FFS Test")
include("FluxFFS.jl")

display("Running Flux FFF Test")
include("FluxFFF.jl")


display("Running Radial Continuum Test")
include("RadialContinuum.jl")

display("Running Radial FSS Test")
include("RadialFSS.jl")

display("Running Radial FFS Test")
include("RadialFFS.jl")

display("Running Radial FFF Test")
include("RadialFFF.jl")
