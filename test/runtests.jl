module MIDTests

#testing seems to take fkn ages, but it will be a good place to write down some expected results.
#probably eventually want to swap this from Axel_q to whatever q-profile we end up using.
#we probbaly want to expand this, but this will confirm that our 3 most basic cases are still giving the same tae as before.
#eg add islands and damping etc.

#in total takes ~3mins to run atm. May want to change to smaller grids for FFF case.
#Probably want all the tests to run at lower grids, more just making sure nothing is broken, we dont actually care about the result

#getting Helmholtz to work is going to be tricky
#think we will need to import the Helmholtz module
#which will need different matrix initialisation etc
#but I think we can multiple dipatch on using real versions of stuff? perhaps?
#probably not actuallly, RIP.
#perhaps we don't worry about Helmholtz, but we just create helmholtz convergence testing with our generic version...


using Test
using MID

#think we should be cleverer about this.
#damping covers non-idea fss and shift and ivert
#in island covers ffs.
#so we can just take an fff case?
#don't think we need the like 10 different version of everthing tbh!

display("Testing continuum calculation")
@time @testset "Continuum" begin include("Continuum.jl") end

#test fff
#and spectrum slicing
display("Testing basic finite elements")
@time @testset "Basic" begin include("Basic.jl") end


#tests island coordinates
#and ffs
display("Testing Island Coordinates")
@time @testset "Island Coordinates" begin include("Island.jl") end

#test fss
#radial coordinate
#shift and invert
display("Testing Continuum Damping")
@time @testset "Damping" begin include("Damping.jl") end

#display("Running QFM Tests")
#@time @testset "QFM" begin include("QFM/QFM.jl") end


#display("Running IO Tests")
#@time @testset "IO" begin include("IO/IO.jl") end



end
