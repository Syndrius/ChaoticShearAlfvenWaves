"""
Base class that just imports everything. We will want a description of the package here and the author etc

# Modules that are still cooked
 - ExtraSpectra -> this will depend on the final state of our verification.
 - IslandContinuum -> Will depend on ExtraSpectra, as we probbaly want to combine into a continuum module.
 - Misc -> Mainly needs a new name and a cleanup, but is a bit subject to change.
 - Plotting -> Pretty cooked, some functions actually need to be fixed, mainly labels are the problem.
 - WeakForm -> Bit cooked, mainly needs a good clean!

 - Fix up plotting
 - Fix up/finish docstrings. Overall code still needs some cleaning, final results will depend on verification via CKA.
 - Should probably check the regularization condition for ϕ
 - Check gadi results and see if we can confirm some convergence
 - Would be good to do some kind of maybe bar graph with contributions from different poloidal mode numbers, would show the island coupling quite nicely I think.
 - determine the width of the islands, much trickier with newest form, either need to approximate with r=0.5 or maybe make a root finding function?? Either way it is pretty cooked. Version on mathematic does not seem to match mode structure, but mode structure seems unreliable for width, should compare with poincare plot.
 - see if we can engineer a case with an island without a TAE so that we can look at the island modes. This may be difficult, as it could be difficult to compare island frequency vs normal continuum frequency.
 - Is it possible to test the theory that island damping is small because tae only interacts for limited r values? I.e frequnecy upshift is across fatest part of island (i think...) which may only be ~10% of the possible values? Ideally it would only be ~1% which might explain our significantly lower damping rate? Perhaps it would be ~10% of θ values and ~10% of ζ values, resulting in a total of ~1%?? Not sure how to argue this or test it or anything tbh.
 - May need to change our profile/setup, seems like any island at all causes overlap, as tae freq is only just above lower gap. -> could be good to have a 3/2 or even 2/1 island or something.
 - We should flip the sign of n, either in general or in island case. having a mix is no bloody good.
 - See if we can find the top of up-shift as a frequency in our global code. i.e to a big island, and see if we can find a frequency that matches what the island continuum code predicts is the top. This does not match very well, hopefully due to not enough mode resolution. -> matches better with Gadi results!
 - combine individual parts of TAE to see if we can see resonance with island, not sure how best to do this... will have to modify our surface function thing! -> i.e. see if we can understand what our eigenfunctions actually mean?
 - Maybe make doctrings look better -> use @doc macro
 - maybe split args and kwargs in doctrings.
 - Integrate the 2d fem stuff properly, and test the damn stuff, -> will want parallel stuf me thinks.
 - Look into profview again, might need to understand views/make use of them. -> may not be as simple as we think! Sounds like sometimes a copy is better, ie if we are using a view of some huge structure for only a small part of it, then the copy can actually be better as it doesn't have to carry everything else with it, or access memory in a weird order.
 - Perhaps make our parallel code more specific to the memory used by each proc, i.e get each proc to work over the indicies it owns? Should be possible with the index_to_grid thingo right? At least to make them mostly on their own memory shouldn't be that bad, will be difficult to deal with the overlap though!
 - So multiple dipatch doesn't work with kwargs, could be time to bin kwargs then! Will be v annoying for solve etc, probably worth keeping them their!
 - Two methods give extremly different results... not good! Probably need to get fem2d working in parallel for proper testing, but we can try the non-damping results, and see if we can replicate the tae freq.
 - may need to double check periodicity on simpler cases
 - maybe we should change problem to accept strings so we can explain the possible options when somehting doesn't work!
 - FFS seems much better at predicting bowden singular, will need to check with a higher resolution though! and also actually check the damping! -> with diagonal met it and lowish res, it seems to be basically perfect! (without damping!) Predicts 0.3258 vs 0.3259, and there is one clear tae (maybe 2!) vs FSS case which has like 10 tae like things. Still obvs need to predict damping!
 - Alot of plotting needs to be fixed, in particular, ms plot, continuum plot have wrong labels.
 - Check Axel's case with FFS, as it is almost perf for bowden singular.
 - Axel's code takes fkn ages to solve, (obvs) but was able to save the matrices, so we can probably solve them in julia??? could be a good idea???
 - To and from file are now cooked, as the grids will depend on the type of solver used, probably need to add a label first, ie GridT: FSS etc
 - Why is parallel mode structure of tae so much better??? Doesn't really make sense...

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



include("Misc/Misc.jl")

using MID.Misc; export matrix_dim



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

