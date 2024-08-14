

"""
    full_spectrum_solve(; Wmat, Imat, grids::GridsT, efuncs=true::Bool, reconstruct=true::Bool, resistivity=false::Bool, R0::Float64)

Uses Julias inbuilt eigen function from LinearAlgebra. Solves the generalised eigenvalue problem Wϕ = ω^2Iϕ, finding the full spectrum of eigenvalues.

### Args
Wmat - W matrix.
Imat - I matrix.
grids::GridsT - Grids used in the problem, used for reconstruction.
efuncs::Bool - Return eigenfunctions with values.
reconstruct::Bool - Whether to reconstruct the eigenfunctions into 3d.
resistivity::Bool - Whether we are solving with resistivity, if not matricies are Hermitian.
R0::Float64 - Major radius, used for normalising.
"""
function full_spectrum_solve(; Wmat, Imat, resistivity=false::Bool)

    if resistivity
        vals, funcs = eigen(Matrix(Wmat), Matrix(Imat))
        return evals, efuncs
    else
        evals, efuncs = eigen(Hermitian(Matrix(Wmat)), Hermitian(Matrix(Imat)))
        return evals, efuncs
    end


end



"""
    arpack_solve(; Wmat, Imat, grids::GridsT, efuncs=true::Bool, nev=20::Int64, σ=0.0::Float64, reconstruct=false::Bool, R0::Float64)

Uses shift and invert to solve for the nev nearest eigenvalues to σ using Arpack. Solves the generalised eigenvalue problem Wϕ = ω^2Iϕ.
Note that Arpack version 0.5.4 is broken but is the default version installed with packagemanger. Install version 0.5.3 with add Arpack@0.5.3

### Args
Wmat - W matrix.
Imat - I matrix.
grids::GridsT - Grids used in the problem, used for reconstruction.
efuncs::Bool - Return eigenfunctions with values.
σ::Float64=0.0 - Find nev nearest evals to σ.
reconstruct::Bool - Whether to reconstruct the eigenfunctions into 3d.
nev::Int64 - Number of eigenvalues to solve for.
R0::Float64 - Major radius, used for normalising.
"""
function arpack_solve(; Wmat, Imat, nev=100::Int64, σ=0.0::Float64, geo::GeoParamsT)

    #un-normalise the target frequency for the shift and invert
    tae_freq = σ^2 / geo.R0^2

    evals, efuncs = eigs(Wmat, Imat, nev=nev, ritzvec=true, sigma=tae_freq)

    return evals, efuncs

end