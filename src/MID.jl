"""
Base class that just imports everything. We will want a description of the package here and the author etc

# Modules that are still cooked
 - Continuum -> Mainly needs to be adjusted to new method of Io. -> and new grid structure -> probably fix once we get island metric invented. -> completly cooked now! -> may need to test this with the island case.


 ###### More Urgent
 - Start using Gadi's version of Julia, currently ours is using fkn heaps of memory
- Map the island results to normal space to see whats going on.
- See if there is any evidence of the gae interaction with the new system
- test m=0 boundaries with fss. Will be annoying.



- Perhaps https://github.com/fredrikekre/Literate.jl if we ever want to share this garbage.
- Probably is not time to change from (r, θ, ζ) into (x1, x2, x3), and change the grids process
think it is best to have a single init_grid function, that just has args for periodicity, boundaries fourier vs fem etc. Alternatively, coords could be labeled as (r, p, t) as in radial poloidal and toroidal? Perhaps that would be clearest?
- pick a lane with met or metric. naming convention of our functions are completly random! And make sure use of ! for in place functions is consistent
- Split spectrum up into construct and solve modules.
- Perhaps it would be nicer to have structure definition of a given module in the main file? just for easier access? Although we are probably supposed to be using the ? help or whatever.
- Start export/using only the functions etc that are actually needed, probably be a lot better. I think ideally this will stop each module from importing the entire module, rather it will import only a few required functions. -> also can make the notation a bit nice, look at QFM for example.
- Perhaps it is time to move plotting to MIDviz?
 - maybe we should change problem to accept strings - or symbols like plot :green etc so we can explain the possible options when something doesn't work!
 - Modify the grids structure, it is horrible, in particular the cont grids is terrible.
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
 - May need to look and see if there is any corelation for the imaginary modes coming from the islands.
 - It could be worth trying ffs, as we only ever seem to need 0, 2, 4, and the negs for ζ. Maybe more efficient???
 - Eventually want to delete all these random af julia files inside the MID folder.
 - Derivative stuff is no longer working.
 - They way we do normalisation is inconsistent and kind of weird...




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

using MID.Geometry; export toroidal_metric!
using MID.Geometry; export no_delta_metric!
using MID.Geometry; export slab_metric!
using MID.Geometry; export island_metric!
using MID.Geometry; export diagonal_toroidal_metric!
using MID.Geometry; export flux_toroidal_metric!
using MID.Geometry; export cylindrical_metric!
using MID.Geometry; export IslandT
using MID.Geometry; export init_island

#using MID.Geometry; export ContIslandT



include("MagneticField/MagneticField.jl")

using MID.MagneticField; export axel_dens
using MID.MagneticField; export Axel_q
using MID.MagneticField; export island_damping_q
using MID.MagneticField; export island_q
using MID.MagneticField; export uniform_dens
using MID.MagneticField; export bowden_singular_dens
using MID.MagneticField; export default_island_q
using MID.MagneticField; export bowden_singular_q
using MID.MagneticField; export comparison_bowden_dens
using MID.MagneticField; export comparison_bowden_q
using MID.MagneticField; export fu_dam_q
using MID.MagneticField; export island_3_2_q
using MID.MagneticField; export symmetric_q
using MID.MagneticField; export flr_q
using MID.MagneticField; export test_q
using MID.MagneticField; export Axel_island_q
using MID.MagneticField; export island_mode_q
using MID.MagneticField; export island_mode_21
using MID.MagneticField; export gae_isl_q
using MID.MagneticField; export gae_isl_dens
using MID.MagneticField; export tae_isl_damping_q
using MID.MagneticField; export chaos_q


include("Structures/Structures.jl")

#seems like we don't have to export the structure, and instead can just use the constructors.
#there is to many things here.
using MID.Structures; export GeoParamsT
using MID.Structures; export ProblemT
using MID.Structures; export GridsT
using MID.Structures; export FLRT
using MID.Structures; export init_grids
using MID.Structures; export init_problem
using MID.Structures; export inst_grids
using MID.Structures; export find_ind
using MID.Structures; export rfem_grid
using MID.Structures; export afem_grid
using MID.Structures; export asm_grid



include("Indexing/Indexing.jl")



include("Basis/Basis.jl")



include("Integration/Integration.jl")


include("PostProcessing/PostProcessing.jl")



include("Io/Io.jl")

using MID.Io; export inputs_to_file
using MID.Io; export inputs_from_file
using MID.Io; export evals_from_file
using MID.Io; export efunc_from_file
using MID.Io; export fortran_process


include("QFM/QFM.jl")

using MID.QFM; export qfm_continuum
using ..QFM; export construct_surfaces
using ..QFM; export farey_tree



include("Plotting/Plotting.jl") 

using MID.Plotting; export potential_plot
using MID.Plotting; export continuum_plot
using MID.Plotting; export contour_plot
using MID.Plotting; export contour_zeta_plot
using MID.Plotting; export surface_plot
using MID.Plotting; export plot_surfs



include("WeakForm/WeakForm.jl")







include("Spectrum/Spectrum.jl")

using MID.Spectrum; export construct
using MID.Spectrum; export arpack_solve
using MID.Spectrum; export full_spectrum_solve
using MID.Spectrum; export compute_spectrum
using MID.Spectrum; export compute_spectrum_qfm
using MID.Spectrum; export spectrum_from_file


#this is now cooked!

include("Continuum/Continuum.jl")

using MID.Continuum; export island_continuum
using MID.Continuum; export island_width
using MID.Continuum; export continuum
using MID.Continuum; export cyl_cont



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

