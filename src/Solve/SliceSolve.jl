"""
    solve(Pmat::SparseMatrixCSC, Qmat::SparseMatrixCSC, solver::SliceSolverT)

Uses the shift and invert method for multiple targets to get multiple 'slices' of the spectrum, allowing a larger portion to be found.
Note that Arpack version 0.5.4 is broken but is the default version installed with packagemanger. Install version 0.5.3 with add Arpack@0.5.3

### Args
- Pmat::SparseMatrixCSC - P matrix.
- Qmat::SparseMatrixCSC - Q matrix.
- sovler::SliceSolverT - solver type that informs the targets for each slice and the number of eigenvalues to get.
"""
function solve(Pmat::SparseMatrixCSC, Qmat::SparseMatrixCSC, solver::SliceSolverT)

    eval_list =  ComplexF64[]
    efunc_list = zeros(ComplexF64, Pmat.m, length(solver.targets)*solver.nev)

    for (i, target) in enumerate(solver.targets)

        evals, efuncs = eigs(Pmat, Qmat, nev=solver.nev, ritzvec=true, sigma=target)

        eval_list = vcat(eval_list, evals)
        efunc_list[:, (i-1)*solver.nev+1:i*solver.nev] = efuncs[:, :]
    end

    return eval_list, efunc_list

end
