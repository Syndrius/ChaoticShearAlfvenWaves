"""
Base class that just imports everything. We will want a description of the package here and the author etc

####################
Current term is back on, it is working fine, may need higher res for best results compared to case without current term.
####################



- Based on fft stuff in qfm, we may need to double check fft in the spectral method to make sure the scaling is actually correct
- Changing to x1,x2,x3 is better long term, but we probably want some functions to still be written in terms of r, θ, ζ, or even ψ, θ, ζ as unique labels
- Add interpolations in properly
- still some room for optimization by removing a few allocations, i.e. in local_to_global.
- Due for a big clean
- Mapping needs hek loads of work.
- QFM Could probably be described properly now, could do with a final fix up.
- Change anal to analytical for proper release
- Maybe remove double ups for sllice solving, and remove boundary evals, not a huge problem typically. -> done for parallel, not for serial...
- Extrapolate the continuum by using Interpolation, shown that it could work, unsure how important this actually is.
- Option of creating grids just by giving an array directly. -> will require significant changes to file storage! 
- We may want to change global basis to two functions, one that just scales the H functions and one that adds the phase factor, as this would be more convenient for Mapping
- Ideally get the continuum functions to return an evalsT struct, may not be possible


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

#things to actually do!
#fix perN continuum
#fix test cases, including helmholtz, island and damping etc -> helmholtz will be no more -> use generic FEM (change name to MID lol) for benchmarking in thesis.
#delete all the extra random af files and stuff.
#get the extra packages actually working
#will be essential for the examples...
#remove gaussquadrature from construct, perhaps fft as well.
#will need to fix q-profiles and density etc -> rename fu_dam_q to quadratic_q
#we could probably at least change some of the indexing to Cartesian indices, eg the test/trial in fff is gross
#get qfm working again, will be fookin annoying. -> might be worth fully fixing up construct etc first, so functions looks the same. -> may even want to split the QFM up a bit.
#islands, the structs and all related ufnctions, are still cooked af.

####################################

#ideally, we can still use the Helmholtz case for testing! -> not going to be practical tbh!
#we might have to settle for the 3d version though!
#but it should still be possible!

include("Structures/Structures.jl")

using ..Structures; export find_ind


#move MetT to structures
#create Geometry struct 
#init will be in here.
include("Geometry/Geometry.jl")

using ..Geometry; export init_geometry 



#change to fields
#move islands to here -> maybe not tbh as islands are needed for the island met...
#move BFieldT to structures.
#create initialisation structure.
include("Fields/Fields.jl")

using ..Fields; export init_fields
using ..Fields; export fu_dam_q
using ..Fields; export qfm_q
using ..Fields; export low_shear_qfm_q
using ..Fields; export qfm_benchmark_q
using ..Fields; export island_q
using ..Fields; export cantori_q



#think grids/basis/integration still kind of go together
#same as geometry/fields/weakform.
#then I guess we have construct/qfm/mapping?

#needs to be more clearly defined as a single module.
include("Basis/Basis.jl")

#merged with grids.
#new version of Gridding is defs better.
include("Grids/Grids.jl")

using ..Grids; export init_grids, init_grid, init_fem_grid, init_sm_grid


#maybe move some fft stuff here? probably not worth it tbh.
#maybe we could shift more of the integration into here
#just to reduce the huge number of loops needed?
include("Integration/Integration.jl")


#almost certianly cooked af.
include("QFM/QFM.jl")

using ..QFM; export construct_surfaces
using ..QFM; export farey_tree
using ..QFM; export lowest_rationals
using ..QFM; export surface_guess
using ..QFM; export compute_jac #perhaps shouldn't be exported



#can probbaly remove types to make this more indep.
include("Io/Io.jl")

using ..Io; export inputs_to_file
using ..Io; export inputs_from_file
using ..Io; export evals_from_file
using ..Io; export efunc_from_file
#using MID.Io; export fortran_process #not sure what to do with this tbh, clearly belongs somewhere else. Will have to see how much it is used in the future


#probably mostly ok.
include("PostProcessing/PostProcessing.jl")


#mostly ok
#think this form is actually better tbh.
#this can go much earlier in the order
#we will also move the problem into the weakform
#as the problem is the key piece that allows the weakform to function!
#ideallly, we should be able to make this work without needing feilds or geometry at all.
#but they will be needed for initialisation.
include("WeakForm/WeakForm.jl")

using ..WeakForm; export init_problem, inst_problem

#need to fix local_to_global and maybe change integration?
include("Construct/Construct.jl")


#should still change from solve to spectrum I think.
include("Solve/Solve.jl")

#using ..Solve; export compute_spectrum
#using ..Solve; export compute_spectrum_qfm
using ..Solve; export init_solver

include("Spectrum/Spectrum.jl")

using ..Spectrum; export compute_spectrum
using ..Spectrum; export analytical_spectrum


#this is a cooked af file
#maybe this can be simplified a bit from other things.
include("Mapping/Mapping.jl") #hopefully not a mistake



end

