"""
    solve(Wmat::SparseMatrixCSC, Imat::SparseMatrixCSC, solver::FullSpectrumSolverT)

Uses Julias inbuilt eigen function from LinearAlgebra. Solves the generalised eigenvalue problem Wϕ = ω^2Iϕ, finding the full spectrum of eigenvalues.

### Args
- Wmat::SparseMatrixCSC - W matrix.
- Imat::SparseMatrixCSC - I matrix.
- solver::FullSpectrumSolverT - solver struct, informs if problem is ideal or not, meaning the matrices are Hermitian or not.
"""
function solve(Wmat::SparseMatrixCSC, Imat::SparseMatrixCSC, solver::FullSpectrumSolverT)

    if solver.ideal
        evals, efuncs = eigen(Hermitian(Matrix(Wmat)), Hermitian(Matrix(Imat)))
    else
        evals, efuncs = eigen(Matrix(Wmat), Matrix(Imat))
    end

    return evals, efuncs

end
