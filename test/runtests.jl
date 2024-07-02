module MIDTests

#testing seems to take fkn ages, but it will be a good place to write down some expected results.
#probably eventually want to swap this from Axel_q to whatever q-profile we end up using.
#we probbaly want to expand this, but this will confirm that our 3 most basic cases are still giving the same tae as before.
#eg add islands and damping etc.

#in total takes ~3mins to run atm. May want to change to smaller grids for FFF case.

using Test
using MID

@time @testset "FSS" begin include("FSSTests/FSSTest.jl") end

@time @testset "FFS" begin include("FFSTests/FFSTest.jl") end

#this takes 2 and half minutes for fk sake, may need to reduce the size significantly.
@time @testset "FFF" begin include("FFFTests/FFFTest.jl") end


end