"""

Main module for finding the 3d spectrum of the TAE wave equation. Includes two main files,
 - Construct, where the two matrices W and I are constructed.
 - Solve, where the generalised eigenvalue problem Wϕ = ω^2 Iϕ is solved.

"""
module Spectrum


using FFTW
using FastGaussQuadrature
using SparseArrays
using LinearAlgebra
using Arpack 
using Printf
#using DelimitedFiles -> don't think this module needs this anymore!


using MID.Geometry
using MID.Basis
using MID.Indexing
using MID.Integration
using MID.MagneticField
using MID.WeakForm
using MID.Inputs


export construct_and_solve
export spectrum_from_file


include("Construct.jl")

export construct


include("Solve.jl")

export arpack_solve
export full_spectrum_solve



"""
    construct_and_solve(; prob::ProblemT, grids::GridsT, efuncs=true::Bool, σ=0.0::Float64, reconstruct=true::Bool, full_spectrum=false::Bool, nev=20::Int64)

Constructs the two matrices and solves. Can either solve the full spectrum via inbuilt solving (slow), or use a shift and invert to find nev amount of eigenvalues closest to σ (fast).

### Args
- prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
- grids::GridT - Grids to solve over.
- efuncs::Bool=true - Return eigenfunctions with values.
- σ::Float64=0.0 - Find nev nearest evals to σ when solving with arpack.
- reconstruct::Bool=true - Whether to reconstruct the eigenfunctions into 3d, needed for plotting. Do not reconstruct if writing efuncs to file.
- full_spectrum::Bool=false - Whether to solve for the full spectrum with inbuilt solver (slow) or use arpack (fast).
- nev::Int64=20 - Number of eigenvalues to solve for if using Arpack.
"""
function construct_and_solve(; prob::ProblemT, grids::GridsT, efuncs=true::Bool, σ=0.0::Float64, reconstruct=true::Bool, full_spectrum=false::Bool, nev=20::Int64)

    display("Constructing...")
    W, I = construct(prob, grids)
    display("Solving...")
    if full_spectrum 
        if prob.δ == 0.0
            ω, ϕ = full_spectrum_solve(Wmat=W, Imat=I, grids=grids, efuncs=efuncs, reconstruct=reconstruct, resistivity=false, R0=prob.geo.R0)
        else
            ω, ϕ = full_spectrum_solve(Wmat=W, Imat=I, grids=grids, efuncs=efuncs, reconstruct=reconstruct, resistivity=true, R0=prob.geo.R0)
        end
    else
        ω, ϕ = arpack_solve(Wmat=W, Imat=I, grids=grids, efuncs=efuncs, nev=nev, σ=σ, reconstruct=reconstruct, R0=prob.geo.R0)
    end
end


"""
    spectrum_from_file(; dir::String, σ::Float64, nev=20::Int64)

Computes the spectrum from inputs read from file, assumes Arpack solve.

# Args
dir::String - Directory where inputs are stored and results will be written to.
σ::Float64 - Target frequency. 
nev=20::Int64 - Number of eigenvalues to solve for.
"""
function spectrum_from_file(; dir::String, σ::Float64, nev=20::Int64)

    prob, grids = inputs_from_file(dir=dir)

    ω, ϕ = construct_and_solve(prob=prob, grids=grids, σ=σ, nev=nev, reconstruct=false)

    eigvals_to_file(ω=ω, filename=dir * "eigvals.dat")
    eigfuncs_to_file(ϕ=ϕ, filename=dir * "eigfuncs.dat")
end

end