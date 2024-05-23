

"""
Uses Julias inbuild eigen function from LinearAlgebra. Solves the generalised eigenvalue problem Wϕ = ω^2Iϕ, finding the full spectrum of eigenvalues.

# Args
- Wmat W matrix.
- Imat I matrix.
- grids::GridsT Grids used in the problem, used for reconstruction.
- efuncs::Bool Return eigenfunctions with values.
- reconstruct::Bool Whether to reconstruct the eigenfunctions into 3d.
- resistivity::Bool Whether we are solving with resistivity, if not matricies are Hermitian.
- R0::Float64 Major radius, used for normalising.
"""
function full_spectrum_solve(;Wmat, Imat, grids::GridsT, efuncs=true::Bool, reconstruct=true::Bool, resistivity=false::Bool, R0::Float64)

    if efuncs
        if resistivity
            vals, funcs = eigen(Matrix(Wmat), Matrix(Imat))
            if reconstruct
                phi = reconstruct_phi(funcs, length(vals), nr, nθ, nζ)
                return R0 .* sqrt.(vals), phi
            else
                return R0 .* sqrt.(vals), funcs
            end
        else
            #may need a reconstruct flag in the future or a better way to write phi to files.
            vals, funcs = eigen(Hermitian(Matrix(Wmat)), Hermitian(Matrix(Imat)))
            if reconstruct
                phi = reconstruct_phi(funcs, length(vals), grids.rd.N, grids.pmd.count, grids.tmd.count)
                return R0 .* sqrt.(vals), phi
            else
                return R0 .* sqrt.(vals), funcs
            end
            
        end
    else
        if resistivity 
            return R0 .* sqrt.(eigvals(Matrix(Wmat), Matrix(Imat)))
        else
            return R0 .* sqrt.(eigvals(Hermitian(Matrix(Wmat)), Hermitian(Matrix(Imat))))
        end
    end


end



"""
Uses shift and invert to solve for the nev nearest eigenvalues to σ using Arpack. Solves the generalised eigenvalue problem Wϕ = ω^2Iϕ.
Note that Arpack version 0.5.4 is broken but is the default version installed with packagemanger. Install version 0.5.3 with add Arpack@0.5.3

# Args
- Wmat W matrix.
- Imat I matrix.
- grids::GridsT Grids used in the problem, used for reconstruction.
- efuncs::Bool Return eigenfunctions with values.
- σ::Float64=0.0 Find nev nearest evals to σ.
- reconstruct::Bool Whether to reconstruct the eigenfunctions into 3d.
- nev::Int64 Number of eigenvalues to solve for.
- R0::Float64 Major radius, used for normalising.
"""
function arpack_solve(; Wmat, Imat, grids::GridsT, efuncs=true::Bool, nev=20::Int64, σ=0.0::Float64, reconstruct=false::Bool, R0::Float64)

    if efuncs
        #which=:LM here seems to cook it in the same way 0.5.4 does
        #wonder why this ever worked???
        vals, funcs = eigs(Wmat, Imat, nev=nev, ritzvec=true, sigma=σ)
        if reconstruct
            phi = reconstruct_phi(funcs, length(vals), grids.rd.N, grids.pmd.count, grids.tmd.count)
            return R0 .* sqrt.(vals), phi
        else
            return R0 .* sqrt.(vals), funcs
        end
    end

    return R0 .* sqrt.(eigs(Wmat, Imat, nev=nev, ritzvec=false, sigma=σ))

end