"""
#TODO
This module is a complete mess and will almost certainly need to change. I think we will just create a continuum module, unsure where convergence tests will go, they may not be as useful anymore
as most convergence tests will need parallel.

Module for auxillary spectra functions we are interested in. In particular contains function for
 - Computing the continuum without islands
 - Computing the spectrum from the two-mode TAE equation derived by Berk et al 1992, for comparison with literature
 - Running convergence tests for damping resutls
 - Running our code with some approximations in built to closer replicate the analytical equations of Berk et al 1992 for benchmarking.

"""

module ExtraSpectra

using MID.Geometry
using MID.Indexing
using MID.MagneticField
using MID.WeakForm
using MID.Structures
using MID.Spectrum


using FFTW
using LinearAlgebra
using FastGaussQuadrature
using Arpack
using SparseArrays
using Printf
using Accessors 
using DelimitedFiles


include("Continuum.jl")

export continuum


include("TwoMode.jl")

export two_mode


include("Convergence.jl")

export convergence_test
export read_convergence_data
export two_mode_convergence

include("Analytical.jl")

export analytical_construct_and_solve

end