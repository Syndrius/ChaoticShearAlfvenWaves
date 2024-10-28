"""
Base class that just imports everything. We will want a description of the package here and the author etc

# Modules that are still cooked
 - Continuum -> Mainly needs to be adjusted to new method of Io. -> and new grid structure -> probably fix once we get island metric invented. -> completly cooked now!
 - WeakForm -> Bit cooked, mainly needs a good clean!
 - Spectrum -> Mainly just construct, bit of a nightmare!


 ###### More Urgent
 - Start using Gadi's version of Julia, currently ours is using fkn heaps of memory
 - Perhaps change the solver to do a range rather than a target, will need to get this to actually work.
 -  Investigate memory usage, Axel used 800x250x10...
 - If we can solve with more evals and less memory usage, will need to paralise post processing.
 - Try the weakform based on Zhisongs equation for the laplacian term.



 - maybe we should change problem to accept strings - or symbols like plot :green etc so we can explain the possible options when something doesn't work!
 - Add try catch to sqrt in solve, most of the time it is just because of ~0 numbers, but it owuld be good to have a warning rather than just always take abs. -> Maybe in this case we dont return the normalised ones? Or should we always have a normalise flag??? -> cka does this better, below some tolerance they are just set to zero.
 - boundary conditions may need modification for flr, and perhaps the m=1 stuff is still not working properly.
 - fix all the examples and extra spectra garbage
 - Reconstruct phi can probbaly be made more efficient -> not sure if we still want to do it one at a time...
 - remove extra exports from this file.
 - Use of kwargs is inconsistent and sometimes annoying.
 - Maybe we should start removing the modes far from the centre as they tend to be garbage.
 - should move to gadi's version of Julia I think. -> we are also running out of space in home directory...
 - Think ideally we would only instantiate grids once...
 - Combine our interpolation thing into MID. (or perhaps make a new package???)
 - Fix Hermite interpolation.
 - Perhaps we shouls always return the periodified version of the potential as this is better for plotting and perhaps interpolation, could just be part of postprocessing
 - May need to look and see if there is any corelation for the imaginary modes coming from the islands.
 - It could be worth trying ffs, as we only ever seem to need 0, 2, 4, and the negs for ζ. Maybe more efficient???
 - Eventually want to delete all these random af julia files inside the MID folder.
 - Derivative stuff is no longer working.
 - WeakForm is seriously cooked, needs to be actually fixed not just cleaned, i.e. stop declaring matrices at every step...
 - They way we do normalisation is inconsistent and kind of weird...
 - We probably need to embark on the great memory purge!




#Long Term fixes
 - Do the interpolation and derivatives etc for all cases.
 - Clean up the q-profiles.
 - Fix continuum -> probably when island met is introduced. Probably remove cont grids.
 - Make user friendly inputs for grids and problems.
 - Change the names for the grid creation functions
 - Change q-profiles to just accept a polynomial coefficeints, that way we don't need a billion, -> provided our island q thing works, this shouldn't be to bad.
 
 

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



include("Io/Io.jl")

using MID.Io; export inputs_to_file
using MID.Io; export inputs_from_file
using MID.Io; export evals_from_file
using MID.Io; export efunc_from_file



include("Plotting/Plotting.jl") 

using MID.Plotting; export potential_plot
using MID.Plotting; export continuum_plot
using MID.Plotting; export contour_plot
using MID.Plotting; export surface_plot



include("WeakForm/WeakForm.jl")



include("PostProcessing/PostProcessing.jl")



include("Spectrum/Spectrum.jl")

using MID.Spectrum; export construct
using MID.Spectrum; export arpack_solve
using MID.Spectrum; export full_spectrum_solve
using MID.Spectrum; export compute_spectrum
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
using MID.Mapping; export mapped_continuum 


#include("ExtraSpectra/ExtraSpectra.jl")

#using MID.ExtraSpectra; export continuum
#much of this is probably garbage!
#using MID.ExtraSpectra; export two_mode
#using MID.ExtraSpectra; export convergence_test
#using MID.ExtraSpectra; export read_convergence_data
#using MID.ExtraSpectra; export two_mode_convergence
#using MID.ExtraSpectra; export analytical_construct_and_solve


end

