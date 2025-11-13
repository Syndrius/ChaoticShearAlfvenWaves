module MIDTests

#testing seems to take fkn ages, but it will be a good place to write down some expected results.
#probably eventually want to swap this from Axel_q to whatever q-profile we end up using.
#we probbaly want to expand this, but this will confirm that our 3 most basic cases are still giving the same tae as before.
#eg add islands and damping etc.

#in total takes ~3mins to run atm. May want to change to smaller grids for FFF case.
#Probably want all the tests to run at lower grids, more just making sure nothing is broken, we dont actually care about the result


using Test
using MID

display("Running Basic Tests")
@time @testset "Basic" begin include("Basic/Basic.jl") end

#display("Running QFM Tests")
#@time @testset "QFM" begin include("QFM/QFM.jl") end

#TODO -> need to be able to find a basic case that works reliably for a test
#display("Running Island Coordinate Tests")
#@time @testset "Island Coordinates" begin include("Island/Island.jl") end

#display("Running IO Tests")
#@time @testset "IO" begin include("IO/IO.jl") end

#TODO
#Damping


end
