module ChaoticShearAlfvenWavesTests

#this is currently probably adequate tbh

using Test
using ChaoticShearAlfvenWaves

display("Running Basis Tests")
@time @testset "Basic" begin include("Basic/Basic.jl") end

display("Running QFM Tests")
@time @testset "QFM" begin include("QFM/QFM.jl") end

display("Running Island Tests")
@time @testset "Island" begin include("Island/Island.jl") end

display("Running IO Tests")
@time @testset "IO" begin include("Io/Io.jl") end

display("Running FLR Tests")
@time @testset "FLR" begin include("FLR/FLR.jl") end




end
