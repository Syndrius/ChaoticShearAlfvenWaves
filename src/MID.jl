"""
Base class that just imports everything. We will want a description of the package here and the author etc

# Modules that are still cooked
 - Continuum -> Mainly needs to be adjusted to new method of Io. -> and new grid structure -> probably fix once we get island metric invented. -> completly cooked now! -> may need to test this with the island case.


 ###### More Urgent
 - Start using Gadi's version of Julia, currently ours is using fkn heaps of memory
- Map the island results to normal space to see whats going on.
- Test island coord case with ffs, fss doesn't seem to give gap modes?? but looks much better, perhaps ffs will be best of both worlds.
- May need to start removing the raw data as we are using stupid amounts of data...



- Perhaps https://github.com/fredrikekre/Literate.jl if we ever want to share this garbage.
- Maybe we should have surfs as part of prob? dealing with parallel surfs is going to be a bit of a disaster! Although having surfs inside prob is a lot of bloat! Perhaps just a string that stores the location? As we only need to create the interpolant, don't actually need the surfaces themselves? Think this could be ok. seems kind of silly for everything else to be stored as is though!
- May want to change the way the q-profile is done. May be better to precompute the values, and just pass in q and dq and floats so that there is no uncertainty, will need to profile.
- Probably want to start removing some things, like og island continuum etc from the module. Bit of a trim tbh., also like fortran processs etc. This shot should just be somewhere else. We can just create a MIDtesting folder or something that deals with all this extra shite.
- One day we want to only have a few densities and q-profiles, for function that are commonly used, perhaps we can have an additional file with a heap of q-profiles etc so that they can be called if needed, but won't make the final cut. Unsure how this will go with parallel tbh. Shouldn't be too hard to just import the q-profile functions tbh. Once we have some examples more solidly written then I think this will be a good idea!
- Probably is not time to change from (r, θ, ζ) into (x1, x2, x3), and change the grids process
think it is best to have a single init_grid function, that just has args for periodicity, boundaries fourier vs fem etc. Alternatively, coords could be labeled as (r, p, t) as in radial poloidal and toroidal? Perhaps that would be clearest?
- pick a lane with met or metric. naming convention of our functions are completly random! And make sure use of ! for in place functions is consistent
- Split spectrum up into construct and solve modules.
- Perhaps it would be nicer to have structure definition of a given module in the main file? just for easier access? Although we are probably supposed to be using the ? help or whatever.
- Start export/using only the functions etc that are actually needed, probably be a lot better. I think ideally this will stop each module from importing the entire module, rather it will import only a few required functions. -> also can make the notation a bit nice, look at QFM for example.
- Perhaps it is time to move plotting to MIDviz?
 - maybe we should change problem to accept strings - or symbols like plot :green etc so we can explain the possible options when something doesn't work!
 - Add try catch to sqrt in solve, most of the time it is just because of ~0 numbers, but it owuld be good to have a warning rather than just always take abs. -> Maybe in this case we dont return the normalised ones? Or should we always have a normalise flag??? -> cka does this better, below some tolerance they are just set to zero.
 - boundary conditions may need modification for flr, and perhaps the m=1 stuff is still not working properly.
 - fix all the examples and extra spectra garbage
 - remove extra exports from this file.
 - Use of kwargs is inconsistent and sometimes annoying.
 - Maybe we should start removing the modes far from the centre as they tend to be garbage.
 - Think ideally we would only instantiate grids once...
 - Combine our interpolation thing into MID. (or perhaps make a new package???)
 - Fix Hermite interpolation.
 - Perhaps we shouls always return the periodified version of the potential as this is better for plotting and perhaps interpolation, could just be part of postprocessing
 - Eventually want to delete all these random af julia files inside the MID folder.
 - Derivative stuff is no longer working.
 - They way we do normalisation is inconsistent and kind of weird...
 - Re do all the example cases, that will give us an idea of how our structures should be organised
 - Ideally, the only exports in this file will be like constructors, less of the actual functions and structs and shite.
 - Unsure if compute spectrum here makes sense, but it seems like it is essentially the main goal of this code? Perhaps it could fo in Solve?




#Long Term fixes
 - Do the interpolation and derivatives etc for all cases.
 - Clean up the q-profiles.
 - Fix continuum -> probably when island met is introduced. Probably remove cont grids.
 - Make user friendly inputs for grids and problems.
 - Change the names for the grid creation functions
 - Change q-profiles to just accept a polynomial coefficeints, that way we don't need a billion, -> provided our island q thing works, this shouldn't be to bad.
 - Probably not worth it, but we could plausibly use some kind of cartesian index for the grid point, then it may be possible to just have as single construct functions. -> Still have to figure out to do fft, and how to 'integrate' for spectral method. But they are getting closer and closer to homogeneous.
 
 

"""


module MID


include("Geometry/Geometry.jl")

#using MID.Geometry; export toroidal_metric!
#using MID.Geometry; export no_delta_metric!
#using MID.Geometry; export slab_metric!
#using MID.Geometry; export island_metric!
#using MID.Geometry; export diagonal_toroidal_metric!
#using MID.Geometry; export flux_toroidal_metric!
#using MID.Geometry; export cylindrical_metric!
#using MID.Geometry; export IslandT
using ..Geometry; export init_island

#using MID.Geometry; export ContIslandT



include("Equilibrium/Equilibrium.jl")

using ..Equilibrium; export fu_dam_q
using ..Equilibrium; export qfm_benchmark_q
#using ..Equilibrium; export Axel_q
#using MID.MagneticField; export axel_dens
#using MID.MagneticField; export island_damping_q
#using MID.MagneticField; export island_q
#using MID.MagneticField; export uniform_dens
#using MID.MagneticField; export bowden_singular_dens
#using MID.MagneticField; export default_island_q
#using MID.MagneticField; export bowden_singular_q
#using MID.MagneticField; export comparison_bowden_dens
#using MID.MagneticField; export comparison_bowden_q
#using MID.MagneticField; export island_3_2_q
#using MID.MagneticField; export symmetric_q
#using MID.MagneticField; export flr_q
#using MID.MagneticField; export test_q
#using MID.MagneticField; export Axel_island_q
#using MID.MagneticField; export island_mode_q
#using MID.MagneticField; export island_mode_21
#using MID.MagneticField; export gae_isl_q
#using MID.MagneticField; export gae_isl_dens
#using MID.MagneticField; export tae_isl_damping_q
#using MID.MagneticField; export chaos_q


include("Structures/Structures.jl")

#seems like we don't have to export the structure, and instead can just use the constructors.
#there is to many things here.
#using MID.Structures; export GeoParamsT
#using MID.Structures; export ProblemT
#using MID.Structures; export GridsT
#using MID.Structures; export FLRT
using ..Structures; export init_grids
using ..Structures; export init_grid
using ..Structures; export init_geo
using ..Structures; export init_flr
using ..Structures; export init_problem
#using MID.Structures; export inst_grids
using ..Structures; export find_ind
#using MID.Structures; export rfem_grid
#using MID.Structures; export afem_grid
#using MID.Structures; export asm_grid



include("Indexing/Indexing.jl")



include("Basis/Basis.jl")



include("Integration/Integration.jl")


include("PostProcessing/PostProcessing.jl")



include("Io/Io.jl")

using ..Io; export inputs_to_file
using ..Io; export inputs_from_file
using ..Io; export evals_from_file
using ..Io; export efunc_from_file
#using MID.Io; export fortran_process


include("QFM/QFM.jl")

#using MID.QFM; export qfm_continuum
using ..QFM; export construct_surfaces
using ..QFM; export farey_tree



#include("Plotting/Plotting.jl") 

#using MID.Plotting; export potential_plot
#using MID.Plotting; export continuum_plot
#using MID.Plotting; export contour_plot
#using MID.Plotting; export contour_zeta_plot
#using MID.Plotting; export surface_plot
#using MID.Plotting; export plot_surfs



include("WeakForm/WeakForm.jl")




include("Construct/Construct.jl")

using ..Construct; export construct


include("Solve/Solve.jl")

using ..Solve; export compute_spectrum
using ..Solve; export compute_spectrum_qfm
using ..Solve; export compute_continuum


#include("Spectrum/Spectrum.jl")

#using MID.Spectrum; export construct
#using MID.Spectrum; export arpack_solve
#using MID.Spectrum; export full_spectrum_solve
#using MID.Spectrum; export compute_spectrum
#using MID.Spectrum; export compute_spectrum_qfm
#using MID.Spectrum; export spectrum_from_file


#this is now cooked!

#needs to be removed from the module and put somewhere else!
include("Continuum/Continuum.jl")

using MID.Continuum; export island_continuum
using MID.Continuum; export island_width
using MID.Continuum; export continuum
using MID.Continuum; export cyl_cont



#perhaps this should be a different module
#just to keep the bloat down!
include("Mapping/Mapping.jl")

using MID.Mapping; export tor_to_isl
using MID.Mapping; export isl_to_tor
using MID.Mapping; export isl_to_tor_continuum 
using MID.Mapping; export tor_to_isl_continuum 


#include("ExtraSpectra/ExtraSpectra.jl")

#using MID.ExtraSpectra; export continuum
#much of this is probably garbage!
#using MID.ExtraSpectra; export two_mode
#using MID.ExtraSpectra; export convergence_test
#using MID.ExtraSpectra; export read_convergence_data
#using MID.ExtraSpectra; export two_mode_convergence
#using MID.ExtraSpectra; export analytical_construct_and_solve


end

