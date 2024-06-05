"""
Base class that just imports everything. We will want a description of the package here and the author etc



 - May want to fix up W.jl... Should be more confident with this now, Worth doing properly so we can be confident about Bs≠0 -> note convergence tests where done with new weakforms, but we should be pretty confident that the new weak form is identical to the old!
 - Fix up plotting
 - Change rdata to store the clustered grid sep points rather than the entire grid. Done but there will be places that this still needs to be fixed. i.e. MIDParallel
 - Fix up/finish docstrings.
 - Island continuum? Not sure if we should do Zhisongs version, or wait for our approximate form to come, once this is done we should be able to finish/fix up numerical methods/MID section of thesis. -> Should be able to re-use Zhisongs form just replace hte parameters and note that it will only be accurate near the island due to taylor expansion.
 - Maybe function to read ith eigenfunction? May become very important with large parallel cases.
 - Probably just need to implment axels case directly...
 - Try damping with g_12 etc = 0, and perhaps the cylindrical metric for the current and I terms? - This essentially does not work, may have to more specifically try and replicate the equation specifically, meaning we will need to derive it fk me.
 - Should probably check the regularization condition for ϕ
 - Compute the continuum, 
 - run some tests with Gadi
 - Start writing the paper... :(

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




#TODO include("IslandContinuum/IslandContinuum.jl")

end

