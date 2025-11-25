"""
    solve(Pmat::SparseMatrixCSC, Qmat::SparseMatrixCSC, solver::ShiftInvertSolverT)

Uses shift and invert to solve for the nev nearest eigenvalues to target frequency using Arpack. Solves the generalised eigenvalue problem PΦ = ω²QΦ.
Note that Arpack version 0.5.4 is broken but is the default version installed with packagemanger. Install version 0.5.3 with add Arpack@0.5.3

### Args
- Pmat::SparseMatrixCSC - P matrix.
- Qmat::SparseMatrixCSC - Q matrix.
- solver::ShitfInvertSolver - solver type, informs the number of eigenvalues to solve for and the target frequency.
"""
function solve(Pmat::SparseMatrixCSC, Qmat::SparseMatrixCSC, solver::ShiftInvertSolverT)

    evals, efuncs = eigs(Pmat, Qmat, nev=solver.nev, ritzvec=true, sigma=solver.target)

    return evals, efuncs

end
