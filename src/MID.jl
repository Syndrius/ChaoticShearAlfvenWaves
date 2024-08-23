"""
Base class that just imports everything. We will want a description of the package here and the author etc

# Modules that are still cooked
 - Continuum -> Mainly needs to be adjusted to new method of Io.
 - WeakForm -> Bit cooked, mainly needs a good clean!
 - Spectrum -> Mainly just construct, bit of a nightmare!


 - We should flip the sign of n, either in general or in island case. having a mix is no bloody good. -> have done this now, I am sure issues will arise lol.
 - pf is cooked again, have to actually modify the mode structure when it is used. Not sure how. (maybe best to ignore.) Seems like the eigenvalues are fine though lol.
 - maybe we should change problem to accept strings - or symbols like plot :green etc so we can explain the possible options when something doesn't work!
 - Add try catch to sqrt in solve, most of the time it is just because of ~0 numbers, but it owuld be good to have a warning rather than just always take abs. -> Maybe in this case we dont return the normalised ones? Or should we always have a normalise flag???
 - boundary conditions may need modification for flr, and perhaps the m=1 stuff is still not working properly.
 - Maybe change count to N...
 - fix all the examples and extra spectra garbage
 - Reconstruct phi can probbaly be made more efficient
 - remove extra exports from this file.
 - Use of kwargs is inconsistent and sometimes annoying.
 - Naming of the plotting functions could be more consistent.
 

"""


module MID


include("Geometry/Geometry.jl")

using MID.Geometry; export toroidal_metric!
using MID.Geometry; export no_delta_metric!
using MID.Geometry; export slab_metric!
using MID.Geometry; export diagonal_toroidal_metric!
using MID.Geometry; export IslandT
using MID.Geometry; export ContIslandT



include("MagneticField/MagneticField.jl")

using MID.MagneticField; export axel_dens
using MID.MagneticField; export Axel_q
using MID.MagneticField; export island_damping_q
using MID.MagneticField; export island_q
using MID.MagneticField; export uniform_dens
using MID.MagneticField; export bowden_singular_dens
using MID.MagneticField; export axel_dens
using MID.MagneticField; export default_island_q
using MID.MagneticField; export bowden_singular_q
using MID.MagneticField; export comparison_bowden_dens
using MID.MagneticField; export comparison_bowden_q
using MID.MagneticField; export fu_dam_q
using MID.MagneticField; export island_3_2_q
using MID.MagneticField; export symmetric_q
using MID.MagneticField; export flr_q



include("Structures/Structures.jl")

#seems like we don't have to export the structure, and instead can just use the constructors.
#there is to many things here.
using MID.Structures; export GeoParamsT
using MID.Structures; export ProblemT
using MID.Structures; export GridsT
using MID.Structures; export FLRT
using MID.Structures; export init_fem_grid
using MID.Structures; export init_sm_grid
using MID.Structures; export init_grids
using MID.Structures; export init_problem
using MID.Structures; export inputs_to_file
using MID.Structures; export inputs_from_file
using MID.Structures; export eigvals_to_file
using MID.Structures; export eigfuncs_to_file
using MID.Structures; export instantiate_grids
using MID.Structures; export find_ind



include("Indexing/Indexing.jl")

using MID.Indexing; export matrix_dim


include("Basis/Basis.jl")



include("Integration/Integration.jl")



include("Io/Io.jl")


using MID.Io; export inputs_to_file
using MID.Io; export inputs_from_file
using MID.Io; export evals_from_file
using MID.Io; export efunc_from_file



include("Plotting/Plotting.jl") 

using MID.Plotting; export plot_potential
using MID.Plotting; export plot_continuum
using MID.Plotting; export contour_plot
using MID.Plotting; export surface_plot



include("WeakForm/WeakForm.jl")



include("Spectrum/Spectrum.jl")

using MID.Spectrum; export construct
using MID.Spectrum; export arpack_solve
using MID.Spectrum; export full_spectrum_solve
using MID.Spectrum; export compute_spectrum
using MID.Spectrum; export post_process
using MID.Spectrum; export spectrum_from_file



include("Continuum/Continuum.jl")

using MID.Continuum; export island_continuum
using MID.Continuum; export island_width
using MID.Continuum; export continuum


#include("ExtraSpectra/ExtraSpectra.jl")

#using MID.ExtraSpectra; export continuum
#much of this is probably garbage!
#using MID.ExtraSpectra; export two_mode
#using MID.ExtraSpectra; export convergence_test
#using MID.ExtraSpectra; export read_convergence_data
#using MID.ExtraSpectra; export two_mode_convergence
#using MID.ExtraSpectra; export analytical_construct_and_solve


end

