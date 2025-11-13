"""
    solve(Wmat::SparseMatrixCSC, Imat::SparseMatrixCSC, solver::ShiftInvertSolverT)

Uses shift and invert to solve for the nev nearest eigenvalues to target frequency using Arpack. Solves the generalised eigenvalue problem Wϕ = ω²Iϕ.
Note that Arpack version 0.5.4 is broken but is the default version installed with packagemanger. Install version 0.5.3 with add Arpack@0.5.3

### Args
- Wmat::SparseMatrixCSC - W matrix.
- Imat::SparseMatrixCSC - I matrix.
- solver::ShitfInvertSolver - solver type, informs the number of eigenvalues to solve for and the target frequency.
"""
function solve(Wmat::SparseMatrixCSC, Imat::SparseMatrixCSC, solver::ShiftInvertSolverT)

    evals, efuncs = eigs(Wmat, Imat, nev=solver.nev, ritzvec=true, sigma=solver.target)

    return evals, efuncs

end
