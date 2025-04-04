"""
Base class that just imports everything. We will want a description of the package here and the author etc

####################
Note that Current term has been switched off! Seems to improve the QFM results. Probably need to investigate!!!
Current term also seems to cook basic toroidal cases. Think there is a mistake somewhere. No current term results look significantly better.
Note that the tests will fail without the current term lol.
####################


 ###### More Urgent
- Map the island results to normal space to see whats going on.
- May need to start removing the raw data as we are using stupid amounts of data...


 - Something is wrong with the current term, (I think!) we need double check the weak form, perhaps using the autodiff/Mathematica to verify the computed values.
 - Add try catch to sqrt in solve, most of the time it is just because of ~0 numbers, but it owuld be good to have a warning rather than just always take abs. -> Maybe in this case we dont return the normalised ones? Or should we always have a normalise flag??? -> cka does this better, below some tolerance they are just set to zero.
 - QFM Module is still a complete disaster.
 - As is weakform, Think we need to try and fix Tj and Tl once and for all. -> could be worth actually making sure the numerical approximation of each term is being computed properly!
 - Remainig modueles are pretty good.




#Long Term fixes
 - Perhaps https://github.com/fredrikekre/Literate.jl if we ever want to share this garbage.
 - Do the interpolation and derivatives etc for all cases. And/or move mapping to a separate module.
 - Re do all the example cases, include proper benchmarks and tests.
 - Use of kwargs is inconsistent and sometimes annoying.
 - May want to change the way the q-profile is done. May be better to precompute the values, and just pass in q and dq and floats so that there is no uncertainty, will need to profile.
 - Maybe we should start removing the modes far from the centre as they tend to be garbage.
 - Maybe we should have surfs as part of prob? dealing with parallel surfs is going to be a bit of a disaster! Although having surfs inside prob is a lot of bloat! Perhaps just a string that stores the location? As we only need to create the interpolant, don't actually need the surfaces themselves? Think this could be ok. seems kind of silly for everything else to be stored as is though!
 - Clean up the q-profiles and density-profiles. We could just have a couple, and move the rest into the profiles folder.
 - Change q-profiles to just accept a polynomial coefficeints, that way we don't need a billion, -> provided our island q thing works, this shouldn't be to bad.
 - Maybe we should change from (r, θ, ζ) into (x1, x2, x3), and change the grids process
 
 

"""


module MID


include("Geometry/Geometry.jl")

using ..Geometry; export init_island



include("Equilibrium/Equilibrium.jl")

using ..Equilibrium; export fu_dam_q
using ..Equilibrium; export qfm_benchmark_q



include("Structures/Structures.jl")

using ..Structures; export init_grids
using ..Structures; export init_grid
using ..Structures; export init_geo
using ..Structures; export init_flr
using ..Structures; export init_problem
using ..Structures; export init_solver
using ..Structures; export find_ind



include("Indexing/Indexing.jl")



include("Basis/Basis.jl")



include("Integration/Integration.jl")



include("PostProcessing/PostProcessing.jl")



include("Io/Io.jl")

using ..Io; export inputs_to_file
using ..Io; export inputs_from_file
using ..Io; export evals_from_file
using ..Io; export efunc_from_file
#using MID.Io; export fortran_process #not sure what to do with this tbh, clearly belongs somewhere else. Will have to see how much it is used in the future



include("QFM/QFM.jl")

using ..QFM; export construct_surfaces
using ..QFM; export farey_tree



include("WeakForm/WeakForm.jl")



include("Construct/Construct.jl")



include("Solve/Solve.jl")

using ..Solve; export compute_spectrum
using ..Solve; export compute_spectrum_qfm
using ..Solve; export compute_continuum



end

