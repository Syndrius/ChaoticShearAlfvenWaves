"""
Base class that just imports everything. We will want a description of the package here and the author etc

####################
Current term is back on, it is working fine, may need higher res for best results compared to case without current term.
####################


 ###### More Urgent
- Map the island results to normal space to see whats going on.
- May need to start removing the raw data as we are using stupid amounts of data...
- Add continuum plotting option to ignore eigenmodes far from validity.


- Based on fft stuff in qfm, we may need to double check fft in the spectral method to make sure the scaling is actually correct
- Changing to x1,x2,x3 is better long term, but we probably want some functions to still be written in terms of r, θ, ζ, or even ψ, θ, ζ as unique labels
- Add interpolations in properly
- still some room for optimization by removing a few allocations, i.e. in local_to_global.
- Due for a big clean
- Mapping needs hek loads of work.
- QFM Could probably be described properly now, could do with a final fix up.
- Change anal to analytical for proper release
- Maybe remove double ups for sllice solving, and remove boundary evals, not a huge problem typically.
- Extrapolate the continuum by using Interpolation, shown that it could work, unsure how important this actually is.
- Option of creating grids just by giving an array directly. -> will require significant changes to file storage! 
- We may want to change global basis to two functions, one that just scales the H functions and one that adds the phase factor, as this would be more convenient for Mapping


#Long Term fixes
 - Perhaps https://github.com/fredrikekre/Literate.jl if we ever want to share this garbage.
 - Re do all the example cases, include proper benchmarks and tests. -> tests are better now, but still not ideal!
 - Perhaps add the option of an iota problem instead? would just be a slightly different compute_B.
 - Add try catch to sqrt in solve, most of the time it is just because of ~0 numbers, but it owuld be good to have a warning rather than just always take abs. -> Maybe in this case we dont return the normalised ones? Or should we always have a normalise flag??? -> cka does this better, below some tolerance they are just set to zero.
 - Use of kwargs is inconsistent and sometimes annoying.
 - May want to change the way the q-profile is done. May be better to precompute the values, and just pass in q and dq and floats so that there is no uncertainty, will need to profile.
 - Throughout code, we have have assumed axissymmetry, this is probably a bad idea, axissymmetry should just be implemented via metric, not in B etc. -> mostly fixed now, but worth checking properly -> this is important for qfm, where, for example, dJ[3] can be non-zero, which is often assumed, as it is zero for most cases.
 - May want to add compute_B to problem, our is actually very restrictive based on the form in terms of r etc. Especially noticable with island case, where q profile is in a different location.
 - Maybe we should start removing the modes far from the centre as they tend to be garbage. -> unsure if this should be done automatically, or they should be kept and we just have an option in plotting?
 - Maybe we should have surfs as part of prob? dealing with parallel surfs is going to be a bit of a disaster! Although having surfs inside prob is a lot of bloat! Perhaps just a string that stores the location? As we only need to create the interpolant, don't actually need the surfaces themselves? Think this could be ok. seems kind of silly for everything else to be stored as is though!
 - Clean up the q-profiles and density-profiles. We could just have a couple, and move the rest into the profiles folder.
 - Change q-profiles to just accept a polynomial coefficeints, that way we don't need a billion, -> provided our island q thing works, this shouldn't be to bad.
 - The order of includes in this file has become very precarious, purpose of some modules is perhaps not clear enough
 - we may benefit from more, but smaller modules one day. This has become especially apparent with Mapping, especially given that the main mapping functions are in post-processing, which could be a mistake.
 - Perhaps move Island.jl into its own module
 - Just like we have a mode label, we might want a way to get the index for a given mode number? Rare usage but who knows
 
 

"""


module MID


include("Geometry/Geometry.jl")

using ..Geometry; export init_island



include("Equilibrium/Equilibrium.jl")

using ..Equilibrium; export fu_dam_q
using ..Equilibrium; export qfm_q
using ..Equilibrium; export low_shear_qfm_q
using ..Equilibrium; export qfm_benchmark_q
using ..Equilibrium; export island_q
using ..Equilibrium; export cantori_q



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



include("QFM/QFM.jl")

using ..QFM; export construct_surfaces
using ..QFM; export farey_tree
using ..QFM; export lowest_rationals
using ..QFM; export surface_guess
using ..QFM; export compute_jac #perhaps shouldn't be exported




include("Io/Io.jl")

using ..Io; export inputs_to_file
using ..Io; export inputs_from_file
using ..Io; export evals_from_file
using ..Io; export efunc_from_file
#using MID.Io; export fortran_process #not sure what to do with this tbh, clearly belongs somewhere else. Will have to see how much it is used in the future



include("PostProcessing/PostProcessing.jl")



include("WeakForm/WeakForm.jl")



include("Construct/Construct.jl")



include("Solve/Solve.jl")

using ..Solve; export compute_spectrum
using ..Solve; export compute_spectrum_qfm
using ..Solve; export compute_continuum



include("Mapping/Mapping.jl") #hopefully not a mistake



end

