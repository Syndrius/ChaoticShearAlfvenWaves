"""
Base class that just imports everything. We will want a description of the package here and the author etc



- Add documentation (Documentor.jl") -> has been initiated.
- May need to change the documentation for multiple dispatch methods. -> eg, we may want a generic construct function placeholder, that gives the general description.
- Create Readme
- add title and descirpiotn to this file
- Fix examples
- Fix all the imports and exports
- Delete all the extra garbage.
- change ζ and φ to be consistent with thesis etc.
- make sure all the random af functions we have still actually work.
- stop the Chaotic etc .QFM.function, i.e. actually import the specific functions used? Perhaps this should be done in general?
- Final run through of each module -> make sure each has a proper description.
- perhaps try and remove any more extra packages from CSAW



"""


module ChaoticShearAlfvenWaves


include("Structures/Structures.jl")

using ..Structures; export find_ind, init_flr, init_island


include("Geometry/Geometry.jl")

using ..Geometry; export init_geometry 


include("Fields/Fields.jl")

using ..Fields; export init_fields 
using ..Fields; export quadratic_q, island_q, damping_q, gae_q, cantori_q
using ..Fields; export uniform_dens, damping_dens, gae_dens


include("Basis/Basis.jl")


include("Grids/Grids.jl")

using ..Grids; export init_grids, init_grid, init_fem_grid, init_sm_grid


include("Integration/Integration.jl")


include("QFM/QFM.jl")


include("Io/Io.jl")

using ..Io; export inputs_to_file, surfaces_to_file, surfaces_from_file
using ..Io; export inputs_from_file
using ..Io; export evals_from_file
using ..Io; export efunc_from_file


include("PostProcessing/PostProcessing.jl")

include("WeakForm/WeakForm.jl")

using ..WeakForm; export init_problem, inst_problem

include("Construct/Construct.jl")


include("Solve/Solve.jl")

using ..Solve; export init_solver


include("Spectrum/Spectrum.jl")

using ..Spectrum; export compute_spectrum
using ..Spectrum; export analytical_spectrum

include("Mapping/Mapping.jl") 

using ..Mapping; export qfm_spectrum_to_tor, isl_spectrum_to_tor, qfm_spectrum_to_isl, tor_spectrum_to_isl


end
