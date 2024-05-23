module Spectrum


using FFTW
using FastGaussQuadrature
using SparseArrays
using LinearAlgebra
using Arpack #this package is completly fked.
using Printf
using DelimitedFiles


using MID.Geometry
using MID.Misc
using MID.MagneticField
using MID.WeakForm
using MID.Inputs


export construct_and_solve
export spectrum_from_file
export solve_from_file
export solve_from_file_from_inputs


include("Construct.jl")

export construct


include("Solve.jl")

export arpack_solve
export full_spectrum_solve



"""
Constructs the two matrices and solves. Can either solve the full spectrum via inbuilt solving (slow), or use a shift and invert to find nev amount of eigenvalues closest to σ (fast).

# Args
- prob::ProblemT Struct containing the functions and parameters that define the problem we are solving
- grids::GridT Grids to solve over.
- efuncs::Bool Return eigenfunctions with values.
- σ::Float64=0.0 Find nev nearest evals to σ when solving with arpack.
- reconstruct::Bool Whether to reconstruct the eigenfunctions into 3d.
- full_spectrum::Bool Whether to solve for the full spectrum with inbuilt solver (slow) or use arpack (fast).
- nev::Int64 Number of eigenvalues to solve for if using Arpack.
"""
function construct_and_solve(; prob::ProblemT, grids::GridsT, efuncs=true::Bool, σ=0.0::Float64, reconstruct=true::Bool, full_spectrum=false::Bool, nev=20::Int64)

    W, I = construct(prob=prob, grids=grids)
    
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

#this is only done under the assumption of ocnvergence so will always use Arpack.
#ie this is essentially just for running convergence tests.
#may want a write to file for the normal constrct and solve func.
function spectrum_from_file(; dir::String, freq::Float64, nev=20::Int64)

    prob, grids = inputs_from_file(dir=dir)

    ω, ϕ = construct_and_solve(prob=prob, grids=grids, σ=freq, nev=nev, reconstruct=false)

    eigvals_to_file(ω=ω, filename=dir * "eigvals.dat")
    eigfuncs_to_file(ϕ=ϕ, filename=dir * "eigfuncs.dat")
end

#useful for interacting with the parralle part which is dodge af
function solve_from_file_from_inputs(; dir::String, σ::Float64)

    file_inds = dir * @sprintf("inds.dat")
    file_data = dir * @sprintf("data.dat")

    rows, cols = eachrow(readdlm(file_inds, ',', Int64))
    Wdata, Idata = eachrow(readdlm(file_data, ',', ComplexF64))

    prob, grids = inputs_from_file(dir=dir)

    ω, ϕ = arpack_solve(Wmat=sparse(rows, cols, Wdata), Imat=sparse(rows, cols, Idata), grids=grids, R0=prob.geo.R0, σ=σ)
    eigvals_to_file(ω=ω, filename=dir * "eigvals.dat")
    eigfuncs_to_file(ϕ=ϕ, filename=dir * "eigfuncs.dat")
end

#useful for interacting with the parralle part which is dodge af
function solve_from_file(; dir::String, grids::GridsT, R0::Float64, σ::Float64)

    file_inds = dir * @sprintf("inds.dat")
    file_data = dir * @sprintf("data.dat")

    rows, cols = eachrow(readdlm(file_inds, ',', Int64))
    Wdata, Idata = eachrow(readdlm(file_data, ',', ComplexF64))

    ω, ϕ = arpack_solve(Wmat=sparse(rows, cols, Wdata), Imat=sparse(rows, cols, Idata), grids=grids, R0=R0, σ=σ)
    eigvals_to_file(ω=ω, filename=dir * "eigvals.dat")
    eigfuncs_to_file(ϕ=ϕ, filename=dir * "eigfuncs.dat")
end

end