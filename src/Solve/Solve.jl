
module Solve

using Printf
using SparseArrays
using LinearAlgebra
using Arpack

using ..Structures
using ..Construct
using ..Basis
using ..PostProcessing




include("FullSpectrumSolve.jl")


include("ShiftInvertSolve.jl")



export compute_spectrum
export compute_continuum

"""
    compute_spectrum(; prob::problemt, grids::gridst, σ=0.0::float64, full_spectrum=false::bool, nev=100::int64)

constructs the two matrices, solves for the eigenvalues and eigenfunctions then processes the results. can either solve the full spectrum via inbuilt solving (slow), or use a shift and invert to find nev amount of eigenvalues closest to σ (fast).

### args
- prob::problemt - struct containing the functions and parameters that define the problem we are solving
- grids::gridt - grids to solve over.
- σ::float64=0.0 - find nev nearest evals to σ when solving with arpack.
- full_spectrum::bool=false - whether to solve for the full spectrum with inbuilt solver (slow) or use arpack (fast).
- nev::int64=100 - number of eigenvalues to solve for if using arpack.
"""

function compute_spectrum(; prob::ProblemT, grids::GridsT, target_freq=0.0::Float64, full_spectrum=false::Bool, nev=100::Int64, deriv=false::Bool)

    t1 = time()

    display("Constructing...")
    @allocated W, I = construct(prob, grids)
    mat_size = matrix_size(grids)
    @printf("Construction of %dx%d matrices complete.\n", mat_size, mat_size)
    display("Solving...")
    if full_spectrum 
        #with no non-ideal effects the matrices are hermitian.
        if prob.flr.δ == 0.0 && prob.flr.ρ_i == 0 && prob.flr.δ_e == 0
            @allocated evals, efuncs = full_spectrum_solve(W, I, true)
        #other wise use a non-hermitian solver.
        else
            evals, efuncs = full_spectrum_solve(Wmat=W, Imat=I, ideal=false)
        end
    else
        #un-normalise the target frequency for the shift and invert
        target_freq = target_freq^2 / prob.geo.R0^2 
        evals, efuncs = shift_invert_solve(Wmat=W, Imat=I, nev=nev, target_freq=target_freq)
    end
    @printf("Solving complete, %d eigenvalues found.\n", length(evals))
    display("Post processing...")

    @allocated evals, ϕ, ϕft = post_process(evals, efuncs, grids, prob.geo, deriv)

    display("Finished.")
    @printf("Total time = %f\n", time() - t1)
    return evals, ϕ, ϕft

end



"""
    compute_continuum(prob::ProblemT, grids::ContGridsT, perN=true::Bool)

Computes the shear Aflven continuum for the non-island case. Has the default option to compute the continuum for each n (toroidal mode number) individually.
"""
function compute_continuum(prob::ProblemT, grids::ContGridsT, perN=true::Bool)

    #we can do this because usual cases (i.e. no island) have no toroidal coupling
    #so we can compute the continuum for each n individually
    if perN
        #_, _, _, _, _, nlist, _ = instantiate_grids(grids)
        nlist = mode_list(grids.ζ)

        ωlist = zeros(grids.r.N, grids.θ.N, grids.ζ.N)

        for (i, n) in enumerate(nlist)
            temp_ζgrid = init_grid(type=:as, start=n, N=1)
            temp_grids = init_grids(grids.r, grids.θ, temp_ζgrid)

            ωlist[:, :, i] = continuum(prob, temp_grids)
        end
    else
        #are we ever actually going to want this??
        ωlist = continuum(prob, grids)
        #need to reformat the output.
    end

    return ωlist

end

#
#can do this in the analytical limit, useful for actually having the mode labels.
function analytical_continuum(prob::ProblemT, grids::ContGridsT)


    rgrid = inst_grid(grids.r)
    mlist = mode_list(grids.θ)
    nlist = mode_list(grids.ζ)

    #can make this the same as evalsT for plotting.
    ω = []
    mode_labs = Tuple{Int, Int}[] 
    rms = []

    for r in rgrid


        q, _ = prob.q(r)

        dens = prob.dens(r)

        for m in mlist

            for n in nlist

                push!(mode_labs, (m, n))

                push!(ω, abs(m/q + n) / sqrt(dens)) #already normalised!

                push!(rms, r)

            end
        end

    end

    return EvalsT(ω, rms, mode_labs)

end



end
