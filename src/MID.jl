"""
Base class that just imports everything. We will want a description of the package here and the author etc

# Modules that are still cooked
 - ExtraSpectra -> this will depend on the final state of our verification.
 - IslandContinuum -> Will depend on ExtraSpectra, as we probbaly want to combine into a continuum module.
 - Plotting -> Pretty cooked, some functions actually need to be fixed, mainly labels are the problem.
 - WeakForm -> Bit cooked, mainly needs a good clean!

 - Fix up/finish docstrings. Overall code still needs some cleaning, final results will depend on verification via CKA. -> use @doc macro
 - see if we can engineer a case with an island without a TAE so that we can look at the island modes. This may be difficult, as it could be difficult to compare island frequency vs normal continuum frequency.
 - Is it possible to test the theory that island damping is small because tae only interacts for limited r values? I.e frequnecy upshift is across fatest part of island (i think...) which may only be ~10% of the possible values? Ideally it would only be ~1% which might explain our significantly lower damping rate? Perhaps it would be ~10% of θ values and ~10% of ζ values, resulting in a total of ~1%?? Not sure how to argue this or test it or anything tbh.
 - May need to change our profile/setup, seems like any island at all causes overlap, as tae freq is only just above lower gap. -> could be good to have a 3/2 or even 2/1 island or something.
 - We should flip the sign of n, either in general or in island case. having a mix is no bloody good.
 - See if we can find the top of up-shift as a frequency in our global code. i.e to a big island, and see if we can find a frequency that matches what the island continuum code predicts is the top. This does not match very well, hopefully due to not enough mode resolution. -> matches better with Gadi results! -> need to redo with ffs if that is significantly better.
 - combine individual parts of TAE to see if we can see resonance with island, not sure how best to do this... will have to modify our surface function thing! -> i.e. see if we can understand what our eigenfunctions actually mean?
 - Perhaps make our parallel code more specific to the memory used by each proc, i.e get each proc to work over the indicies it owns? Should be possible with the index_to_grid thingo right? At least to make them mostly on their own memory shouldn't be that bad, will be difficult to deal with the overlap though!
 - Two methods give extremly different results... not good! Probably need to get fem2d working in parallel for proper testing, but we can try the non-damping results, and see if we can replicate the tae freq.
 - maybe we should change problem to accept strings so we can explain the possible options when something doesn't work!
 - FFS seems much better at predicting bowden singular, will need to check with a higher resolution though! and also actually check the damping! -> with diagonal met it and lowish res, it seems to be basically perfect! (without damping!) Predicts 0.3258 vs 0.3259, and there is one clear tae (maybe 2!) vs FSS case which has like 10 tae like things. Still obvs need to predict damping!
 - Alot of plotting needs to be fixed, in particular, ms plot, continuum plot have wrong labels.
 - Check Axel's case with FFS, as it is almost perf for bowden singular.
 - Why is parallel mode structure of tae so much better??? Doesn't really make sense...
 - maybe have a constructing and solving output as well!
 - Probably should add the periodic part to any final plots, I think that will make our periodicity clearer.
 - Add try catch to sqrt in solve, most of the time it is just because of ~0 numbers, but it owuld be good to have a warning rather than just always take abs. -> Maybe in this case we dont return the normalised ones? Or should we always have a normalise flag???
 - May need to change how efuncs are read in, ffs case is already having problemo's. May need to do a single n at a time. Or if we do FFF, we will want a particular slice of ζ I guess.
 - Wonder if we are seeing small damping because n0=4 and tae is n=2?? Think it is time to change our test case unfort. Changing to FFF would probably also fix this!
 - consider normal damping case without island but with mulitple n's. This may highlight that the tae does not get damped by n ≠2 modes even if the continuum is overlapping. May add some evidence that our theory that we need an n=2 island for an n=2 tae to experience significant damping.
 - Similarly, we probbaly need to understand what the different n's mean in the island case.
 - Put all our thoughts onto paper, ideally with pictures of our justification of each. Eg show island influence, show damping not being effected by different n's in no island case etc.
 - Reonstruct phi could probably be sped up by skipping every 2/4/8 indexes.
 - Plot continuum is probably the most cooked function going around!
 - Boundaries are being added twice... not sure why... -> this may be because of the two θbasis funcs, then in 3d it is added 4 times for the two θ and two ζ basis functions. -> has no effect on the frequency. I think this would just scale those basis functions by 2, but because they are zero it doesn't matter. Would be a problemo if we had non-zero boundaries.
 - Maybe create a grid type for continuum, i.e. just the spectral method??
 - Move Gauss quadrature to integration.jl
 - Why is ffs case spiking at r=0 in island cases?


 Two things to try:
    - New q-profile/setup so that n of island matches n of tae, hoepfully this will offer more interaction -> not sure this is possible, unless island is not located at gap, this seems problomatic, as then we cannot use the continuum as verifaction for the damping. May be worth a go though!
    - change to FFF as this seems easier to underttand, then we can consider ζ cross section over different parts of the island.

Two ways to increase the damping/interaction
    - Make the tae and the island have the same toroidal mode number, will require modification, and don't think it is possible to have the island at the gap. -> doesn't seem to work.
    - Make the interaction between the neighbouring toroidal mode (i.e. + n0) and the continuum stronger.


 - Think we are closer to understanding gaps, but we still need to think!
 - Perhaps it would be worth directly implementing a two mode case without any simplifications, see if it matches our larger code. Hopefully not to much effort???



"""


module MID


include("Geometry/Geometry.jl")

using MID.Geometry; export toroidal_metric!
using MID.Geometry; export no_delta_metric!
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



include("Inputs/Inputs.jl")

#seems like we don't have to export the structure, and instead can just use the constructors.
using MID.Inputs; export GeoParamsT
using MID.Inputs; export ProblemT
using MID.Inputs; export GridsT
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
using MID.Plotting; export plot_potential
using MID.Plotting; export plot_sum_potential
using MID.Plotting; export find_ind
using MID.Plotting; export plot_continuum
using MID.Plotting; export plot_phi_surface
using MID.Plotting; export construct_surface
using MID.Plotting; export mode_structure



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



include("ExtraSpectra/ExtraSpectra.jl")

using MID.ExtraSpectra; export continuum
#much of this is probably garbage!
using MID.ExtraSpectra; export two_mode
using MID.ExtraSpectra; export convergence_test
using MID.ExtraSpectra; export read_convergence_data
using MID.ExtraSpectra; export two_mode_convergence
using MID.ExtraSpectra; export analytical_construct_and_solve




include("IslandContinuum/IslandContinuum.jl")

using MID.IslandContinuum; export island_continuum
using MID.IslandContinuum; export island_width

end

