"""
    solve(Pmat::SparseMatrixCSC, Qmat::SparseMatrixCSC, solver::FullSpectrumSolverT)

Uses Julias inbuilt eigen function from LinearAlgebra. Solves the generalised eigenvalue problem PΦ = ω^2QΦ, finding the full spectrum of eigenvalues.

### Args
- Pmat::SparseMatrixCSC - P matrix.
- Qmat::SparseMatrixCSC - Q matrix.
- solver::FullSpectrumSolverT - solver struct, informs if problem is ideal or not, meaning the matrices are Hermitian or not.
"""
function solve(Pmat::SparseMatrixCSC, Qmat::SparseMatrixCSC, solver::FullSpectrumSolverT)

    if solver.ideal
        evals, efuncs = eigen(Hermitian(Matrix(Pmat)), Hermitian(Matrix(Qmat)))
    else
        evals, efuncs = eigen(Matrix(Pmat), Matrix(Qmat))
    end

    return evals, efuncs

end
