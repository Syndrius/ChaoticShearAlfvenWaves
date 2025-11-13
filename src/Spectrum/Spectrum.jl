"""

This modules solve the eigenvalues problem Wϕ - ω²Iϕ given the two matrices W and I. 
Solving can be done
 - directly via Julia's LinearAlgebra, which is slow and not practical for large grids.
 - using shift and invert to target a specific frequency, obtaining the nev::Int64 nearest eigenvalues.
 - using a 'slicing' method where the shift invert method is used multiple times to build up a larger portion of the spectrum.

Also contains the compute_spectrum function, which is essentially the main() of this package. This function constructs the matrices, solves the eigenvalue problem then processes the output.
"""
module Spectrum

using Printf
#using SparseArrays
#using LinearAlgebra
using ..WeakForm #why do we need this..
using ..Grids
using ..Structures
#import ..ChaoticShearAlfvenWaves: GridsT, SolverT, ProblemT #think this isa ctually a bad idea
#here we are specifically solving this particular problem, I think...
#taking this approach cooks the QFM and Cont stuff. Perhaos there are better places for them?
#we might be able to remove weakform with a different QFM problemT but the grids will be tricky.
using ..Solve
using ..Construct
using ..PostProcessing
using ..Io


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
function compute_spectrum(; prob::ProblemT, grids::GridsT, solver::SolverT, deriv=false::Bool)

    t1 = time()

    display("Constructing...")
    @allocated W, I = construct(prob, grids)
    mat_size = matrix_size(grids)
    @printf("Construction of %dx%d matrices complete.\n", mat_size, mat_size)

    display("Solving...")
    evals, efuncs = solve(W, I, solver)
    @printf("Solving complete, %d eigenvalues found.\n", length(evals))

    display("Post processing...")

    evals, ϕ, ϕft = post_process(evals, efuncs, grids, prob.geo, deriv)

    display("Finished.")
    @printf("Total time = %f\n", time() - t1)
    return evals, ϕ, ϕft

end

#=

"""
    compute_spectrum_qfm(; prob::ProblemT, grids::GridsT, solver::SolverT, surfs::Array{QFMSurfaceT})

Same as compute_spectrum, but uses qfm surfaces.

### args
- prob::ProblemT - struct containing the functions and parameters that define the problem we are solving
- grids::GridT - grids to solve over.
- solver::SolverT - struct storing the parameters for solving.
- surfs::Array{QFMSurfaceT} - array of the qfm surfaces.
"""
function compute_spectrum_qfm(; prob::ProblemT, grids::GridsT, solver::SolverT, surfs::Array{QFMSurfaceT}, deriv=false::Bool)

    t1 = time()

    display("Constructing...")
    @allocated W, I = construct(prob, grids, surfs)
    mat_size = matrix_size(grids)
    @printf("Construction of %dx%d matrices complete.\n", mat_size, mat_size)

    display("Solving...")
    evals, efuncs = solve(W, I, solver)
    @printf("Solving complete, %d eigenvalues found.\n", length(evals))

    display("Post Processing...")
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
        nlist = mode_list(grids.x3)

        ωlist = zeros(grids.x1.N, grids.x2.N, grids.x3.N)

        for (i, n) in enumerate(nlist)
            temp_x3grid = init_grid(type=:as, start=n, N=1)
            temp_grids = init_grids(grids.x1, grids.x2, temp_x3grid)

            ωlist[:, :, i] = continuum(prob, temp_grids, surfs)
        end
    else
        #are we ever actually going to want this??
        ωlist = continuum(prob, grids, surfs)
        #need to reformat the output.
    end

    return ωlist

end
=#

"""
    compute_continuum(prob::ProblemT, grids::ContGridsT, perN=true::Bool)

Computes the shear Aflven continuum for the non-island case. Has the default option to compute the continuum for each n (toroidal mode number) individually.
"""
#TODO, changing to evals struct means the perN thing won't work anymore!
function compute_spectrum(prob::ProblemT, grids::ContGridsT, perN=false::Bool)

    #we can do this because usual cases (i.e. no island) have no toroidal coupling
    #so we can compute the continuum for each n individually
    #doesn't currently work!
    if perN
        #_, _, _, _, _, nlist, _ = instantiate_grids(grids)
        nlist = mode_list(grids.x3)

        ωlist = zeros(grids.x1.N, grids.x2.N, grids.x3.N)

        for (i, n) in enumerate(nlist)
            temp_x3grid = init_grid(type=:as, start=n, N=1)
            temp_grids = init_grids(grids.x1, grids.x2, temp_x3grid)

            ωlist[:, :, i] = continuum(prob, temp_grids)
        end
    else
        #are we ever actually going to want this??
        ωlist = continuum(prob, grids)
        #need to reformat the output.
    end

    evals = post_process(ωlist, grids, prob.geo)

    return evals

end


"""
    analytical_continuum(prob::ProblemT, grids::ContGridsT)

Computes the continuum analytically, assuming cylindrical geometry and no perturbations. Useful as a comparison.
"""
function analytical_spectrum(prob::ProblemT, grids::ContGridsT)


    x1grid = inst_grid(grids.x1)
    mlist = mode_list(grids.x2)
    nlist = mode_list(grids.x3)

    ω = []
    mode_labs = Tuple{Int, Int}[] 
    x1ms = []

    for x1 in x1grid
        q, _ = prob.q(x1)

        dens = prob.dens(x1)

        for m in mlist, n in nlist

            push!(mode_labs, (m, n))

            push!(ω, abs(m/q + n) / sqrt(dens)) #already normalised!

            push!(x1ms, x1)
        end

    end

    return EvalsT(ω, x1ms, mode_labs)

end


#basic functionality works here, need more serious testing with actual anon-functions
function spectrum_from_file(dir::String, deriv::Bool=false)
    uninst_prob, grids, solver = inputs_from_file(dir=dir)

    prob = Structures.inst_problem(uninst_prob)

    evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver, deriv=deriv)

    evals_to_file(evals, dir)

    mkpath(dir * "/efuncs")
    mkpath(dir * "/efuncs_ft")

    efuncs_to_file(ϕ, ϕft, dir)

end

end
