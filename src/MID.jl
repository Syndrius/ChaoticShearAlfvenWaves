"""
Base class that just imports everything. We will want a description of the package here and the author etc



 - Fix up plotting
 - Fix up/finish docstrings.
 - Maybe function to read ith eigenfunction? May become very important with large parallel cases.
 - Should probably check the regularization condition for ϕ
 - Check gadi results and see if we can confirm some convergence
 - Make island continuum good, and see if we can figure out what the cut-off island A should be, and see if it matches our code.
 - Would be good to do some kind of maybe bar graph with contributions from different poloidal mode numbers, would show the island coupling quite nicely I think.
 - determine the width of the islands, much trickier with newest form, either need to approximate with r=0.5 or maybe make a root finding function?? Either way it is pretty cooked.
 - see if we can engineer a case with an island without a TAE so that we can look at the island modes. This may be difficult, as it could be difficult to compare island frequency vs normal continuum frequency.
 - Is it possible to test the theory that island damping is small because tae only interacts for limited r values? I.e frequnecy upshift is across fatest part of island (i think...) which may only be ~10% of the possible values? Ideally it would only be ~1% which might explain our significantly lower damping rate? Perhaps it would be ~10% of θ values and ~10% of ζ values, resulting in a total of ~1%?? Not sure how to argue this or test it or anything tbh.
 - May need to change our profile/setup, seems like any island at all causes overlap, as tae freq is only just above lower gap.
 - We should flip the sign of n, either in general or in island case. having a mix is no bloody good.

"""


module MID


include("Geometry/Geometry.jl")

using MID.Geometry; export toroidal_metric!
using MID.Geometry; export no_delta_metric!
using MID.Geometry; export diagonal_toroidal_metric!
using MID.Geometry; export IslandT



include("MagneticField/MagneticField.jl")

using MID.MagneticField; export axel_dens
using MID.MagneticField; export Axel_q
using MID.MagneticField; export island_damping_q
using MID.MagneticField; export island_q
using MID.MagneticField; export uniform_dens
using MID.MagneticField; export bowden_singular_dens
using MID.MagneticField; export axel_dens
using MID.MagneticField; export default_island_q
using MID.MagneticField; export singular_bowden_q
using MID.MagneticField; export comparison_bowden_dens
using MID.MagneticField; export comparison_bowden_q
using MID.MagneticField; export fu_dam_q



include("Misc/Misc.jl")

#most of this has to be exported to allow MIDParallel to work. Would it be better to 
using MID.Misc; export clustered_grid #may not need to share this anymore.
#using MID.Misc; export spectral_grid
#using MID.Misc; export hermite_basis
#using MID.Misc; export compute_boundary_inds
#using MID.Misc; export local_to_global
#using MID.Misc; export create_local_basis




include("Inputs/Inputs.jl")

#seems like we don't have to export the structure, and instead can just use the constructors.
using MID.Inputs; export GeoParamsT
using MID.Inputs; export ProblemT
using MID.Inputs; export GridsT
using MID.Inputs; export init_grids
using MID.Inputs; export init_problem
using MID.Inputs; export inputs_to_file
using MID.Inputs; export inputs_from_file
using MID.Inputs; export eigvals_to_file
using MID.Inputs; export eigfuncs_to_file



include("Plotting/Plotting.jl") #bit of a disaster atm!

using MID.Plotting; export reconstruct_continuum
using MID.Plotting; export plot_potential
using MID.Plotting; export find_ind
using MID.Plotting; export plot_continuum



include("WeakForm/WeakForm.jl")



include("Spectrum/Spectrum.jl")

using MID.Spectrum; export construct
using MID.Spectrum; export arpack_solve
using MID.Spectrum; export full_spectrum_solve
using MID.Spectrum; export construct_and_solve
using MID.Spectrum; export reconstruct_phi
using MID.Spectrum; export spectrum_from_file
using MID.Spectrum; export solve_from_file
using MID.Spectrum; export solve_from_file_from_inputs


include("ExtraSpectra/ExtraSpectra.jl")

using MID.ExtraSpectra; export continuum
using MID.ExtraSpectra; export two_mode
using MID.ExtraSpectra; export convergence_test
using MID.ExtraSpectra; export read_convergence_data
using MID.ExtraSpectra; export two_mode_convergence
using MID.ExtraSpectra; export analytical_construct_and_solve




include("IslandContinuum/IslandContinuum.jl")

using MID.IslandContinuum; export trapped_continuum
using MID.IslandContinuum; export passing_continuum
using MID.IslandContinuum; export island_continuum
using MID.IslandContinuum; export ContIslandT

end

