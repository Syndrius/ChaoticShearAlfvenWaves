"""
Base class that just imports everything. We will want a description of the package here and the author etc


 - get this running on Pc for island convergence tests
 - see if we can verify the code now, may need to compare to original 1d FEM case from Berk. Think Bowdens' singular finite elements code may be a good shout!

 - Need some more examples of code working
 - Need some benchmarking
 - May want to fix up W.jl...
 - Otherwise this package should be kind of complete.
 - Can't say the same for the others tho lol.


# Plots backend is probably cooked
# Arpack version needs to be verified, looks like broken version has been removed from github, so should be fine.
"""


module MID


include("Geometry/Geometry.jl")

using MID.Geometry; export toroidal_metric!
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



include("Misc/Misc.jl")

#most of this has to be exported to allow MIDParallel to work. Would it be better to 
using MID.Misc; export clustered_grid
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



include("Continuum/Continuum.jl")

using MID.Continuum; export continuum



include("Spectrum/Spectrum.jl")

using MID.Spectrum; export construct
using MID.Spectrum; export arpack_solve
using MID.Spectrum; export full_spectrum_solve
using MID.Spectrum; export construct_and_solve
using MID.Spectrum; export reconstruct_phi
using MID.Spectrum; export spectrum_from_file
using MID.Spectrum; export solve_from_file
using MID.Spectrum; export solve_from_file_from_inputs





#solving is too cooked to use atm! May try to integrate later, but way to hard atm.
#gone back to separate package
#include("Parallel/Parallel.jl")

#using MID.Parallel; export par_construct_and_solve
#using MID.Parallel; export par_construct
#using MID.Parallel; export par_construct_to_file

#using MID.Parallel; export par_solve
#using MID.Parallel; export par_spectrum_from_file



#TODO include("IslandContinuum/IslandContinuum.jl")

end

