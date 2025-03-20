"""

Main module for finding the 3d spectrum of the TAE wave equation. Includes three main files,
 - Construct, where the two matrices W and I are constructed.
 - Solve, where the generalised eigenvalue problem Wϕ = ω^2 Iϕ is solved.
 - PostProcess, where the solutions are returned to the physical 3d form, eigenvalues are normalised and the continuum is reconstructed.

"""
module Spectrum


using FFTW
using FastGaussQuadrature
using SparseArrays
using LinearAlgebra
using Arpack 
using Printf



using MID.Geometry
using MID.Basis
using MID.Indexing
using MID.Integration
using MID.MagneticField
using MID.WeakForm
using MID.Structures
using MID.Io
using MID.PostProcessing
using MID.QFM


export compute_spectrum
export compute_spectrum_qfm
export spectrum_from_file


include("Construct.jl")

export construct


include("QFMConstruct.jl")

export construct


include("Solve.jl")

export arpack_solve
export full_spectrum_solve



"""
    compute_spectrum(; prob::ProblemT, grids::GridsT, σ=0.0::Float64, full_spectrum=false::Bool, nev=100::Int64)

Constructs the two matrices, solves for the eigenvalues and eigenfunctions then processes the results. Can either solve the full spectrum via inbuilt solving (slow), or use a shift and invert to find nev amount of eigenvalues closest to σ (fast).

### Args
- prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
- grids::GridT - Grids to solve over.
- σ::Float64=0.0 - Find nev nearest evals to σ when solving with arpack.
- full_spectrum::Bool=false - Whether to solve for the full spectrum with inbuilt solver (slow) or use arpack (fast).
- nev::Int64=100 - Number of eigenvalues to solve for if using Arpack.
"""
function compute_spectrum(; prob::ProblemT, grids::GridsT, target_freq=0.0::Float64, full_spectrum=false::Bool, nev=100::Int64, deriv=false::Bool)

    t1 = time()

    display("Constructing...")
    @allocated W, I = construct(prob, grids)
    mat_size = matrix_size(grids)
    @printf("Construction of %dx%d matrices complete.\n", mat_size, mat_size)
    display("Solving...")
    if full_spectrum 
        #with no non-ideal effects the matrices are Hermitian.
        if prob.flr.δ == 0.0 && prob.flr.ρ_i == 0 && prob.flr.δ_e == 0
            @allocated evals, efuncs = full_spectrum_solve(Wmat=W, Imat=I, ideal=true)
        #other wise use a non-hermitian solver.
        else
            evals, efuncs = full_spectrum_solve(Wmat=W, Imat=I, ideal=false)
        end
    else
        #un-normalise the target frequency for the shift and invert
        target_freq = target_freq^2 / prob.geo.R0^2 
        evals, efuncs = arpack_solve(Wmat=W, Imat=I, nev=nev, target_freq=target_freq)
    end
    @printf("Solving complete, %d eigenvalues found.\n", length(evals))
    display("Post Processing...")

    @allocated evals, ϕ, ϕft = post_process(evals, efuncs, grids, prob.geo, deriv)

    display("Finished.")
    @printf("total time = %f\n", time() - t1)
    return evals, ϕ, ϕft

end


#perhaps we should stop using generic names quite as much, instead we could call this compute spectrum qfm.
#this is causing some precompilation issues, may need to change the name.
function compute_spectrum_qfm(; prob::ProblemT, grids::GridsT, surfs::Array{QFMSurfaceT}, target_freq=0.0::Float64, full_spectrum=false::Bool, nev=100::Int64, deriv=false::Bool)

    t1 = time()

    display("Constructing...")
    @allocated W, I = construct(prob, grids, surfs)
    mat_size = matrix_size(grids)
    @printf("Construction of %dx%d matrices complete.\n", mat_size, mat_size)
    display("Solving...")
    if full_spectrum 
        #with no non-ideal effects the matrices are Hermitian.
        if prob.flr.δ == 0.0 && prob.flr.ρ_i == 0 && prob.flr.δ_e == 0
            @allocated evals, efuncs = full_spectrum_solve(Wmat=W, Imat=I, ideal=true)
        #other wise use a non-hermitian solver.
        else
            evals, efuncs = full_spectrum_solve(Wmat=W, Imat=I, ideal=false)
        end
    else
        #un-normalise the target frequency for the shift and invert
        target_freq = target_freq^2 / prob.geo.R0^2 
        evals, efuncs = arpack_solve(Wmat=W, Imat=I, nev=nev, target_freq=target_freq)
    end
    @printf("Solving complete, %d eigenvalues found.\n", length(evals))
    display("Post Processing...")

    @allocated evals, ϕ, ϕft = post_process(evals, efuncs, grids, prob.geo, deriv)

    display("Finished.")
    @printf("total time = %f\n", time() - t1)
    return evals, ϕ, ϕft

end


"""
    spectrum_from_file(; dir::String, σ::Float64, nev=100::Int64)

Computes the spectrum from inputs read from file, assumes Arpack solve. Writes output to files in same directory.

# Args
dir::String - Directory where inputs are stored and results will be written to.
σ::Float64 - Target frequency. 
nev=100::Int64 - Number of eigenvalues to solve for.
"""
function spectrum_from_file(; dir::String, target_freq::Float64, nev=100::Int64)

    prob, grids = inputs_from_file(dir=dir)

    evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, target_freq=target_freq, nev=nev)

    evals_to_file(evals, dir)
    efuncs_to_file(ϕ, ϕft, dir)
end

end
