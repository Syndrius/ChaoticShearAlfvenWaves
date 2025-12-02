"""
Package for computing the spectrum of Shear Alfven Waves in fusion plasmas with chaotic magentic fields.
This packages solved the generalised eigenvalue problem of the reduced ideal MHD equations for shear Alfven waves.

Written by Matthew Thomas.
See Thomas et al. 2026 for a decription.
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

using ..Io; export inputs_to_file, surfaces_to_file, surfaces_from_file, inputs_from_file, evals_from_file, efunc_from_file


include("PostProcessing/PostProcessing.jl")


include("WeakForm/WeakForm.jl")

using ..WeakForm; export init_problem


include("Construct/Construct.jl")


include("Solve/Solve.jl")

using ..Solve; export init_solver


include("Spectrum/Spectrum.jl")

using ..Spectrum; export compute_spectrum, analytical_spectrum


include("Mapping/Mapping.jl") 

using ..Mapping; export qfm_spectrum_to_tor, isl_spectrum_to_tor, qfm_spectrum_to_isl, tor_spectrum_to_isl


end
