"""
    full_spectrum_solve(; Wmat, Imat, ideal=true::Bool)

Uses Julias inbuilt eigen function from LinearAlgebra. Solves the generalised eigenvalue problem Wϕ = ω^2Iϕ, finding the full spectrum of eigenvalues.

### Args
- Wmat::SparseMatrixCSC - W matrix.
- Imat::SparseMatrixCSC - I matrix.
- ideal::Bool - Whether we are solving with ideal mhd, if so matricies are Hermitian.
"""
function solve(Wmat::SparseMatrixCSC, Imat::SparseMatrixCSC, solver::FullSpectrumSolverT)

    if solver.ideal
        evals, efuncs = eigen(Hermitian(Matrix(Wmat)), Hermitian(Matrix(Imat)))
    else
        evals, efuncs = eigen(Matrix(Wmat), Matrix(Imat))
    end

    return evals, efuncs

end
