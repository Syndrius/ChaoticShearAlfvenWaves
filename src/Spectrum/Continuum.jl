
"""
    compute_continuum(prob::ProblemT, grids::ContGridsT)

Computes the shear Aflven continuum for the non-island case.
"""
function compute_spectrum(prob::ProblemT, grids::ContGridsT)

    ωlist = continuum(prob, grids)

    evals = post_process(ωlist, grids, prob.geo)

    return evals

end


"""
    analytical_continuum(prob::ProblemT, grids::ContGridsT)

Computes the continuum analytically, assuming cylindrical geometry and no perturbations.
"""
function analytical_spectrum(prob::ProblemT, grids::ContGridsT)


    x1grid = inst_grid(grids.x1)
    mlist = mode_list(grids.x2)
    nlist = mode_list(grids.x3)

    ω = []
    mode_labs = Tuple{Int, Int}[] 
    x1ms = []

    for x1 in x1grid
        q, _ = prob.fields.q(x1)

        dens = prob.fields.dens(x1)

        for m in mlist, n in nlist

            push!(mode_labs, (m, n))

            push!(ω, abs(m/q + n) / sqrt(dens)) #already normalised

            push!(x1ms, x1)
        end

    end

    return EvalsT(ω, x1ms, mode_labs)

end

