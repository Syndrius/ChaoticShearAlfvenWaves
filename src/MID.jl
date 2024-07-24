"""
Base class that just imports everything. We will want a description of the package here and the author etc

# Modules that are still cooked
 - ExtraSpectra -> this will depend on the final state of our verification.
 - IslandContinuum -> Will depend on ExtraSpectra, as we probbaly want to combine into a continuum module.
 - Plotting -> Pretty cooked, some functions actually need to be fixed, mainly labels are the problem.
 - WeakForm -> Bit cooked, mainly needs a good clean!

 - Fix up/finish docstrings. Overall code still needs some cleaning, final results will depend on verification via CKA. -> use @doc macro
 - We should flip the sign of n, either in general or in island case. having a mix is no bloody good.
 - Perhaps make our parallel code more specific to the memory used by each proc, i.e get each proc to work over the indicies it owns? Should be possible with the index_to_grid thingo right? At least to make them mostly on their own memory shouldn't be that bad, will be difficult to deal with the overlap though!
 - maybe we should change problem to accept strings - or symbols like plot :green etc so we can explain the possible options when something doesn't work!
 - Alot of plotting needs to be fixed, in particular, ms plot, continuum plot have wrong labels.
 - Check Axel's case with FFS, as it is almost perf for bowden singular.
 - Why is parallel mode structure of tae so much better??? Doesn't really make sense...
 - maybe have a constructing and solving output as well!
 - Probably should add the periodic part to any final plots, I think that will make our periodicity clearer.
 - Add try catch to sqrt in solve, most of the time it is just because of ~0 numbers, but it owuld be good to have a warning rather than just always take abs. -> Maybe in this case we dont return the normalised ones? Or should we always have a normalise flag???
 - May need to change how efuncs are read in, ffs case is already having problemo's. May need to do a single n at a time. Or if we do FFF, we will want a particular slice of ζ I guess.
 - Reonstruct phi could probably be sped up by skipping every 2/4/8 indexes.
 - Plot continuum is probably the most cooked function going around!
 - Boundaries are being added twice... not sure why... -> this may be because of the two θbasis funcs, then in 3d it is added 4 times for the two θ and two ζ basis functions. -> has no effect on the frequency. I think this would just scale those basis functions by 2, but because they are zero it doesn't matter. Would be a problemo if we had non-zero boundaries.
 - see if we can estimate how much memory we should be using, and make sure it is not way more than that! -> would help with gadi inputs as well.
 - Maybe a reconstruct continuum from file, so that we can see the continuum for bigger files!
 - change filename to savefile in all cases for clarity.
 - boundary conditions may need modification for flr, and perhaps the m=1 stuff is still not working properly.

 - Really could do with some even very basic tae identification, ie make sure that a specific n (ie of the tae) is actually the maximum, and perhaps make sure that the major m's of the tae are like at least 50% of the max or something? -> ffs is less of a problem, the tae stays a bit more stable when an island is introduced
 - With new mode structure method, we could have a single reconstruct phi, just need to make a potential_size(grids) function.
 - Looks like we need to implement KAW/finite Epar corrections for this to actually work... fk me.
 - Consider integrals.jl as alternative numerical integration
 - Also consider julia --track-allocation=user i.e. run julia with track allocations to get more detailed data on memory allocations.


 Two things to try:
    - New q-profile/setup so that n of island matches n of tae, hoepfully this will offer more interaction -> not sure this is possible, unless island is not located at gap, this seems problomatic, as then we cannot use the continuum as verifaction for the damping. May be worth a go though!
    - change to FFF as this seems easier to underttand, then we can consider ζ cross section over different parts of the island.

Two ways to increase the damping/interaction
    - Make the tae and the island have the same toroidal mode number, will require modification, and don't think it is possible to have the island at the gap. -> doesn't seem to work.
    - Make the interaction between the neighbouring toroidal mode (i.e. + n0) and the continuum stronger.


 - Think we are closer to understanding gaps, but we still need to think!
 - Perhaps it would be worth directly implementing a two mode case without any simplifications, see if it matches our larger code. Hopefully not to much effort???
 - Need to try and find a gae again I think!



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



include("Inputs/Inputs.jl")

#seems like we don't have to export the structure, and instead can just use the constructors.
using MID.Inputs; export GeoParamsT
using MID.Inputs; export ProblemT
using MID.Inputs; export GridsT
using MID.Inputs; export FLRT
using MID.Inputs; export init_fem_grid
using MID.Inputs; export init_sm_grid
using MID.Inputs; export init_grids
using MID.Inputs; export init_problem
using MID.Inputs; export inputs_to_file
using MID.Inputs; export inputs_from_file
using MID.Inputs; export eigvals_to_file
using MID.Inputs; export eigfuncs_to_file



include("Indexing/Indexing.jl")

using MID.Indexing; export matrix_dim



include("Basis/Basis.jl")

using MID.Basis; export hermite_basis
using MID.Basis; export create_local_basis!



include("Integration/Integration.jl")

using MID.Integration; export gauss_integrate



include("Io/Io.jl")

using MID.Io; export inputs_from_file
using MID.Io; export inputs_to_file
using MID.Io; export eigvals_to_file
using MID.Io; export eigfuncs_to_file



include("Plotting/Plotting.jl") #bit of a disaster atm!

using MID.Plotting; export reconstruct_continuum
using MID.Plotting; export reconstruct_slab_continuum
using MID.Plotting; export plot_contour_poincare
using MID.Plotting; export reconstruct_continuum_n
using MID.Plotting; export plot_potential
using MID.Plotting; export plot_sum_potential
using MID.Plotting; export find_ind
using MID.Plotting; export plot_continuum
using MID.Plotting; export plot_phi_surface
using MID.Plotting; export construct_surface
using MID.Plotting; export mode_structure
using MID.Plotting; export contour_plot



include("WeakForm/WeakForm.jl")



include("Spectrum/Spectrum.jl")

using MID.Spectrum; export construct
using MID.Spectrum; export construct_zf
using MID.Spectrum; export arpack_solve
using MID.Spectrum; export full_spectrum_solve
using MID.Spectrum; export construct_and_solve
using MID.Spectrum; export reconstruct_phi
using MID.Spectrum; export spectrum_from_file
using MID.Spectrum; export solve_from_file
using MID.Spectrum; export solve_from_file_from_inputs



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

