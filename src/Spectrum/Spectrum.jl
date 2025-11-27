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

using ..Grids
using ..Structures
using ..Solve
using ..Construct
using ..WeakForm
using ..PostProcessing
using ..Io

export compute_spectrum
export analytical_spectrum


include("Continuum.jl")

include("QFM.jl")


"""
    compute_spectrum(prob::ProblemT, grids::GridsT, solver::SolverT; deriv=false::Bool)

Constructs the two matrices, solves for the eigenvalues and eigenfunctions then processes the results. 

### args
- prob::ProblemT - struct containing the functions and parameters that define the problem we are solving
- grids::GridT - grids to solve over.
- solver::SolverT - struct storing the parameters for solving.
- deriv::Bool=false - If derivatives of eigenfunctions should be returned.
"""
function compute_spectrum(prob::ProblemT, grids::GridsT, solver::SolverT; deriv=false::Bool)

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


"""
    compute_spectrum(dir::String; deriv=false)

Computes the spectrum from inputs written to file.
"""
function compute_spectrum(dir::String; deriv=false)
    uninst_prob, grids, solver = inputs_from_file(dir)

    #this function is quite awkward because jld2 cannot write anonymous functions to file
    #so they have to be recreated.
    prob = WeakForm.inst_problem(uninst_prob.fields, uninst_prob.geo, uninst_prob.flr)

    evals, ϕ, ϕft = compute_spectrum(prob, grids, solver, deriv=deriv)

    evals_to_file(evals, dir)

    mkpath(dir * "/efuncs")
    mkpath(dir * "/efuncs_ft")

    efuncs_to_file(ϕ, ϕft, dir)

end

end
