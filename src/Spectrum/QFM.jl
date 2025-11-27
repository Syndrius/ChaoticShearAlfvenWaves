
"""
    compute_spectrum(prob::ProblemT, grids::GridsT, solver::SolverT, surfs::Array{QFMSurfaceT}, deriv=false::Bool)

Computes the spectrum using QFM surfaces.

### args
- prob::ProblemT - struct containing the functions and parameters that define the problem we are solving
- grids::GridT - grids to solve over.
- solver::SolverT - struct storing the parameters for solving.
- surfs::Array{QFMSurfaceT} - array of the qfm surfaces.
- deriv::Bool=false - If derivatives of eigenfunctions should be returned.
"""
function compute_spectrum(prob::ProblemT, grids::GridsT, solver::SolverT, surfs::Array{QFMSurfaceT}; deriv=false::Bool)

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
    compute_continuum(prob::ProblemT, grids::ContGridsT, surfs::Array{QFMSurfaceT)

Computes the shear Aflven continuum with QFM surfaces.
"""
function compute_spectrum(prob::ProblemT, grids::ContGridsT, surfs::Array{QFMSurfaceT})

    ωlist = continuum(prob, grids, surfs)

    evals = post_process(ωlist, grids, prob.geo)

    return evals

end


"""
Computes the spectrum from inputs written to file.
"""
function compute_spectrum(dir::String, surfs::Array{QFMSurfaceT}; deriv=false)
    uninst_prob, grids, solver = inputs_from_file(dir)

    prob = WeakForm.inst_problem(uninst_prob.fields, uninst_prob.geo, uninst_prob.flr)

    evals, ϕ, ϕft = compute_spectrum(prob, grids, solver, surfs, deriv=deriv)

    evals_to_file(evals, dir)

    mkpath(dir * "/efuncs")
    mkpath(dir * "/efuncs_ft")

    efuncs_to_file(ϕ, ϕft, dir)

end
