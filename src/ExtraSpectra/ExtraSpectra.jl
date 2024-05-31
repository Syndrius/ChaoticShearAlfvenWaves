"""

Module for auxillary spectra functions we are interested in. In particular contains function for
 - Computing the continuum without islands
 - Computing the spectrum from the two-mode TAE equation derived by Berk et al 1992, for comparison with literature
 - Running convergence tests for damping resutls

"""

module ExtraSpectra

using MID.Geometry
using MID.Misc
using MID.MagneticField
using MID.WeakForm
using MID.Inputs
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


end