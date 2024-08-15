
"""
    full_spectrum_solve(; Wmat, Imat, resistivity=false::Bool)

Uses Julias inbuilt eigen function from LinearAlgebra. Solves the generalised eigenvalue problem Wϕ = ω^2Iϕ, finding the full spectrum of eigenvalues.

### Args
- Wmat::SparseMatrixCSC - W matrix.
- Imat::SparseMatrixCSC - I matrix.
- ideal::Bool - Whether we are solving with ideal mhd, if so matricies are Hermitian.
"""
function full_spectrum_solve(; Wmat::SparseMatrixCSC, Imat::SparseMatrixCSC, ideal=true::Bool)

    if ideal
        evals, efuncs = eigen(Matrix(Wmat), Matrix(Imat))
    else
        evals, efuncs = eigen(Hermitian(Matrix(Wmat)), Hermitian(Matrix(Imat)))
    end

    return evals, efuncs

end



"""
    arpack_solve(; Wmat, Imat, nev=100::Int64, σ=0.0::Float64, geo::GeoParamsT)

Uses shift and invert to solve for the nev nearest eigenvalues to σ using Arpack. Solves the generalised eigenvalue problem Wϕ = ω^2Iϕ.
Note that Arpack version 0.5.4 is broken but is the default version installed with packagemanger. Install version 0.5.3 with add Arpack@0.5.3

### Args
- Wmat::SparseMatrixCSC - W matrix.
- Imat::SparseMatrixCSC - I matrix.
- target_freq::Float64=0.0 - Find nev nearest evals to σ.
- nev::Int64 - Number of eigenvalues to solve for.
- geo::GeoParamsT - Geometry struct, used for un-normalising target frequency.
"""
function arpack_solve(; Wmat::SparseMatrixCSC, Imat::SparseMatrixCSC, nev=100::Int64, target_freq=0.0::Float64, geo::GeoParamsT)

    #un-normalise the target frequency for the shift and invert
    tae_freq = target_freq^2 / geo.R0^2

    evals, efuncs = eigs(Wmat, Imat, nev=nev, ritzvec=true, sigma=tae_freq)

    return evals, efuncs

end