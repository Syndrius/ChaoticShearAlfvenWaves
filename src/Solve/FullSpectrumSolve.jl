"""
    full_spectrum_solve(; Wmat, Imat, ideal=true::Bool)

Uses Julias inbuilt eigen function from LinearAlgebra. Solves the generalised eigenvalue problem Wϕ = ω^2Iϕ, finding the full spectrum of eigenvalues.

### Args
- Wmat::SparseMatrixCSC - W matrix.
- Imat::SparseMatrixCSC - I matrix.
- ideal::Bool - Whether we are solving with ideal mhd, if so matricies are Hermitian.
"""
function full_spectrum_solve(; Wmat::SparseMatrixCSC, Imat::SparseMatrixCSC, ideal=true::Bool)
#function full_spectrum_solve(; Wmat, Imat, ideal=true::Bool)

    if ideal
        evals, efuncs = eigen(Hermitian(Matrix(Wmat)), Hermitian(Matrix(Imat)))
    else
        evals, efuncs = eigen(Matrix(Wmat), Matrix(Imat))
    end

    return evals, efuncs

end
