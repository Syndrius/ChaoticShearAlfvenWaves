"""

This modules solve the eigenvalues problem Wϕ - ω²Iϕ given the two matrices W and I. 
Solving can be done
 - directly via Julia's LinearAlgebra, which is slow and not practical for large grids.
 - using shift and invert to target a specific frequency, obtaining the nev::Int64 nearest eigenvalues.
 - using a 'slicing' method where the shift invert method is used multiple times to build up a larger portion of the spectrum.

Also contains the compute_spectrum function, which is essentially the main() of this package. This function constructs the matrices, solves the eigenvalue problem then processes the output.
"""
module Solve

using Printf
using SparseArrays
using LinearAlgebra
using Arpack

using ..Structures
using ..Construct
using ..Basis
using ..PostProcessing
using ..QFM



include("FullSpectrumSolve.jl")


include("ShiftInvertSolve.jl")


include("SliceSolve.jl")


export compute_spectrum
export compute_spectrum_qfm
export compute_continuum


"""
    compute_spectrum(; prob::ProblemT, grids::GridsT, solver::SolverT)

Constructs the two matrices, solves for the eigenvalues and eigenfunctions then processes the results. 

### args
- prob::ProblemT - struct containing the functions and parameters that define the problem we are solving
- grids::GridT - grids to solve over.
- solver::SolverT - struct storing the parameters for solving.
"""
function compute_spectrum(; prob::ProblemT, grids::GridsT, solver::SolverT)

    t1 = time()

    display("Constructing...")
    @allocated W, I = construct(prob, grids)
    mat_size = matrix_size(grids)
    @printf("Construction of %dx%d matrices complete.\n", mat_size, mat_size)
    display("Solving...")
    #with no non-ideal effects the matrices are hermitian.
    evals, efuncs = solve(W, I, solver)
    @printf("Solving complete, %d eigenvalues found.\n", length(evals))
    display("Post processing...")

    deriv = false # not ideal, not sure if we will ever bother with this again
    evals, ϕ, ϕft = post_process(evals, efuncs, grids, prob.geo, deriv)

    display("Finished.")
    @printf("Total time = %f\n", time() - t1)
    return evals, ϕ, ϕft

end

"""
    compute_spectrum_qfm(; prob::ProblemT, grids::GridsT, solver::SolverT, surfs::Array{QFMSurfaceT})

Same as compute_spectrum, but uses qfm surfaces.

### args
- prob::ProblemT - struct containing the functions and parameters that define the problem we are solving
- grids::GridT - grids to solve over.
- solver::SolverT - struct storing the parameters for solving.
- surfs::Array{QFMSurfaceT} - array of the qfm surfaces.
"""
function compute_spectrum_qfm(; prob::ProblemT, grids::GridsT, solver::SolverT, surfs::Array{QFMSurfaceT})

    t1 = time()

    display("Constructing...")
    @allocated W, I = construct(prob, grids, surfs)
    mat_size = matrix_size(grids)
    @printf("Construction of %dx%d matrices complete.\n", mat_size, mat_size)

    display("Solving...")
    evals, efuncs = solve(W, I, solver)
    @printf("Solving complete, %d eigenvalues found.\n", length(evals))

    display("Post Processing...")
    deriv = false # not ideal, not sure if we will ever bother with this again
    @allocated evals, ϕ, ϕft = post_process(evals, efuncs, grids, prob.geo, deriv)
    display("Finished.")

    @printf("total time = %f\n", time() - t1)
    return evals, ϕ, ϕft

end


"""
    compute_continuum(prob::ProblemT, grids::ContGridsT, perN=true::Bool)

Computes the shear Aflven continuum for the non-island case. Has the default option to compute the continuum for each n (toroidal mode number) individually.
"""
function compute_continuum(prob::ProblemT, grids::ContGridsT, surfs::Array{QFMSurfaceT}, perN=true::Bool)

    #we can do this because usual cases (i.e. no island) have no toroidal coupling
    #so we can compute the continuum for each n individually
    if perN
        #_, _, _, _, _, nlist, _ = instantiate_grids(grids)
        nlist = mode_list(grids.ζ)

        ωlist = zeros(grids.r.N, grids.θ.N, grids.ζ.N)

        for (i, n) in enumerate(nlist)
            temp_ζgrid = init_grid(type=:as, start=n, N=1)
            temp_grids = init_grids(grids.r, grids.θ, temp_ζgrid)

            ωlist[:, :, i] = continuum(prob, temp_grids, surfs)
        end
    else
        #are we ever actually going to want this??
        ωlist = continuum(prob, grids, surfs)
        #need to reformat the output.
    end

    return ωlist

end

"""
    compute_continuum(prob::ProblemT, grids::ContGridsT, perN=true::Bool)

Computes the shear Aflven continuum for the non-island case. Has the default option to compute the continuum for each n (toroidal mode number) individually.
"""
function compute_continuum(prob::ProblemT, grids::ContGridsT, perN=true::Bool)

    #we can do this because usual cases (i.e. no island) have no toroidal coupling
    #so we can compute the continuum for each n individually
    if perN
        #_, _, _, _, _, nlist, _ = instantiate_grids(grids)
        nlist = mode_list(grids.ζ)

        ωlist = zeros(grids.r.N, grids.θ.N, grids.ζ.N)

        for (i, n) in enumerate(nlist)
            temp_ζgrid = init_grid(type=:as, start=n, N=1)
            temp_grids = init_grids(grids.r, grids.θ, temp_ζgrid)

            ωlist[:, :, i] = continuum(prob, temp_grids)
        end
    else
        #are we ever actually going to want this??
        ωlist = continuum(prob, grids)
        #need to reformat the output.
    end

    return ωlist

end


"""
    analytical_continuum(prob::ProblemT, grids::ContGridsT)

Computes the continuum analytically, assuming cylindrical geometry and no perturbations. Useful as a comparison.
"""
function analytical_continuum(prob::ProblemT, grids::ContGridsT)


    rgrid = inst_grid(grids.r)
    mlist = mode_list(grids.θ)
    nlist = mode_list(grids.ζ)

    ω = []
    mode_labs = Tuple{Int, Int}[] 
    rms = []

    for r in rgrid
        q, _ = prob.q(r)

        dens = prob.dens(r)

        for m in mlist, n in nlist

            push!(mode_labs, (m, n))

            push!(ω, abs(m/q + n) / sqrt(dens)) #already normalised!

            push!(rms, r)
        end

    end

    return EvalsT(ω, rms, mode_labs)

end


end
