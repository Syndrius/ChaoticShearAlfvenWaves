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

#submodules to fix
#construct
#mapping


#things to actually do!
#fix perN continuum
#fix test cases, including helmholtz, island and damping etc -> helmholtz will be no more -> use generic FEM (change name to MID lol) for benchmarking in thesis.
#delete all the extra random af files and stuff.
#get the extra packages actually working
#will be essential for the examples...
#get qfm working again, will be fookin annoying. -> might be worth fully fixing up construct etc first, so functions looks the same. -> may even want to split the QFM up a bit.
#islands, the structs and all related ufnctions, are still cooked af.
#looks like we might need to fix slepcwrap after all lol -> cannot really get MIDParallel to compile anymore. old MPI version is conflicting with ordinary Diff eq.
#alternatively, we can just create two different environments, one with parallel and one with ordinarydiffeq, until it is fixed. -> may be a more practical solution until thesis is done.
#allow single island input in init_fields..
#change W, I to P, Q to match paper and eventually thesis.
#perhaps even change ζ to φ to match paper -> ideally we can determine what toroidal and cylindrical will actually be in general and make our code consistent with that!
#gotta change the fkn spelling of separatrix....

####################################

#ideally, we can still use the Helmholtz case for testing! -> not going to be practical tbh!
#we might have to settle for the 3d version though!
#but it should still be possible!

#good, assuming nothing else changes.
include("Structures/Structures.jl")

using ..Structures; export find_ind, init_flr, init_island


#good
include("Geometry/Geometry.jl")

using ..Geometry; export init_geometry 


#good
include("Fields/Fields.jl")

using ..Fields; export init_fields 
using ..Fields; export quadratic_q, island_q, damping_q, gae_q, cantori_q
using ..Fields; export uniform_dens, damping_dens, gae_dens


#good
include("Basis/Basis.jl")


#good
include("Grids/Grids.jl")

using ..Grids; export init_grids, init_grid, init_fem_grid, init_sm_grid


#good
include("Integration/Integration.jl")


#good
include("QFM/QFM.jl")


#good
include("Io/Io.jl")

using ..Io; export inputs_to_file
using ..Io; export inputs_from_file
using ..Io; export evals_from_file
using ..Io; export efunc_from_file
#using MID.Io; export fortran_process #not sure what to do with this tbh, clearly belongs somewhere else. Will have to see how much it is used in the future


#good
include("PostProcessing/PostProcessing.jl")

#good
include("WeakForm/WeakForm.jl")

using ..WeakForm; export init_problem, inst_problem

#need to fix local_to_global and maybe change integration?
include("Construct/Construct.jl")


#good
include("Solve/Solve.jl")

using ..Solve; export init_solver


#good
include("Spectrum/Spectrum.jl")

using ..Spectrum; export compute_spectrum
using ..Spectrum; export analytical_spectrum


#this is a cooked af file
#maybe this can be simplified a bit from other things.
include("Mapping/Mapping.jl") #hopefully not a mistake



end

