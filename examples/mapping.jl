
#example of mapping perhaps, or just include within qfm

dir = abspath(joinpath(pathof(MID), "../../test/data/"))

surfs = load_object(dir * "benchmark_surfaces.jld2");
