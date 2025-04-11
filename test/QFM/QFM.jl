
display("Computing surfaces")

R0=4.0

#define the non-resonant island
k = 0.005
isl = init_island(m0=3, n0=2, A=k)

geo = init_geo(R0=R0)

prob = init_problem(q=fu_dam_q, geo=geo, isl=isl)

qlist, plist = farey_tree(3, 2, 1, 3, 1)

guess_list = 0.5 .* ones(length(qlist));

surfs = construct_surfaces(plist, qlist, guess_list, prob);

display("Running Continuum Test")
include("Continuum.jl")

display("Running FSS Test")
include("FSS.jl")

display("Running FFS Test")
include("FFS.jl")

display("Running FFF Test")
include("FFF.jl")

