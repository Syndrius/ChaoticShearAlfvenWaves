
module PostProcessing

using FFTW
using EllipsisNotation

using MID.Structures
using MID.Indexing
using MID.Basis


export post_process


include("Reconstruction.jl")


include("FT.jl")


include("Continuum.jl")



"""
    post_process(evals::AbstractArray, efuncs::Array{ComplexF64}, grids::GridsT, geo::GeoParamsT, deriv=false::Bool)

Processes the eigenvalues and eigenfunctions into useful forms. This involves recosntructing the 1d eigenfunction into the 3d grid representation for ϕ, and taking the fourier transform in θ and ζ. The fourier transform is used to label each eigenvalue with its dominant mode and radial location. The eigenvalues are also normalised to the Alfven frequency.
"""
function post_process(evals::AbstractArray, efuncs::Array{ComplexF64}, grids::GridsT, geo::GeoParamsT, deriv=false::Bool)

    #deriv is only really valid for fff for now, and probably for a while lol.

    #allocates global arrays for storing everyeigenvalue
    ϕ_g, ϕft_g = allocate_phi_arrays(grids, length(evals), deriv=deriv)
    #allocates placeholder arrays used during computation.
    ϕp, ϕpft = allocate_phi_arrays(grids, deriv=deriv)

    #array to store the radial value of each eigenvalue.
    rms = Array{Float64}(undef, length(evals))

    #array to store the mode labels
    mode_labs = Tuple{Int, Int}[] 

    #create a plan for fourier transforming.
    plan = create_ft_plan(ϕpft, grids)

    rgrid = inst_grids(grids)[1]

    #arrays to store the maximum value of ϕft and the correspondning r value.
    rmarray = Array{Int64}(undef, grids.θ.N, grids.ζ.N)
    ϕmarray = Array{Float64}(undef, grids.θ.N, grids.ζ.N)

    
    for i in 1:length(evals)

        #converts the eigenfunction back into the 3d ϕ, including its fourier transformation in θ and ζ
        reconstruct_phi!(efuncs[:, i], grids, ϕp, ϕpft, plan)

        #finds the dominant mode and its radial location for continuum reconstruction.
        rind, mode_lab = label_mode(ϕpft, grids, rmarray, ϕmarray)

        ϕ_g[i, ..] = ϕp
        ϕft_g[i, ..] = ϕpft

        rms[i] = rgrid[rind]
        push!(mode_labs, mode_lab)
    end

    #should only take abs for ideal cases, this is just to avoid
    #tiny negatives that should be zero
    #alterantively, as per cka, just set any tiny negatives to zero.
    #ω = geo.R0 .* sqrt.(abs.(evals))
    ω = geo.R0 .* sqrt.(evals)

    evals = EvalsT(ω, rms, mode_labs)

    return evals, ϕ_g, ϕft_g

end


"""
    allocate_phi_arrays(grids::FSSGridsT, nevals=0; deriv=false)

Allocates empty arrays used for computing phi and its fourier transformation. If nevals is non-zero creates the global phi structures that store all eigenvalues.
"""
function allocate_phi_arrays(grids::FSSGridsT, nevals=0; deriv=false)

    #extends the domain for ifft cases for smoother plotting.
    θsize = ifft_size(grids.θ)
    ζsize = ifft_size(grids.ζ)

    if deriv

        if nevals == 0
            ϕ = Array{ComplexF64}(undef, grids.r.N, θsize, ζsize, 8)
            ϕft = Array{ComplexF64}(undef, grids.r.N, grids.θ.N, grids.ζ.N, 8)
        else
            ϕ = Array{ComplexF64}(undef, nevals, grids.r.N, θsize, ζsize, 8)
            ϕft = Array{ComplexF64}(undef, nevals, grids.r.N, grids.θ.N, grids.ζ.N, 8)
        end
        #placeholder memory for each eigenfunction.
        
    else
        if nevals == 0
            ϕ = Array{ComplexF64}(undef, grids.r.N, θsize, ζsize)
            ϕft = Array{ComplexF64}(undef, grids.r.N, grids.θ.N, grids.ζ.N)
        else
            ϕ = Array{ComplexF64}(undef, nevals, grids.r.N, θsize, ζsize)
            ϕft = Array{ComplexF64}(undef, nevals, grids.r.N, grids.θ.N, grids.ζ.N)
        end
        
    end

    return ϕ, ϕft

end



"""
    allocate_phi_arrays(grids::FFSGridsT, nevals=0; deriv=false)

Allocates empty arrays used for computing phi and its fourier transformation. If nevals is non-zero creates the global phi structures that store all eigenvalues.
"""
function allocate_phi_arrays(grids::FFSGridsT, nevals=0; deriv=false)

    ζsize = ifft_size(grids.ζ)

    if deriv

        if nevals == 0
            ϕ = Array{ComplexF64}(undef, grids.r.N, grids.θ.N, ζsize, 8)
            ϕft = Array{ComplexF64}(undef, grids.r.N, grids.θ.N, grids.ζ.N, 8)
        else
            ϕ = Array{ComplexF64}(undef, nevals, grids.r.N, grids.θ.N, ζsize, 8)
            ϕft = Array{ComplexF64}(undef, nevals, grids.r.N, grids.θ.N, grids.ζ.N, 8)
        end
        
    else
        if nevals == 0
            ϕ = Array{ComplexF64}(undef, grids.r.N, grids.θ.N, ζsize)
            ϕft = Array{ComplexF64}(undef, grids.r.N, grids.θ.N, grids.ζ.N)
        else
            ϕ = Array{ComplexF64}(undef, nevals, grids.r.N, grids.θ.N, ζsize)
            ϕft = Array{ComplexF64}(undef, nevals, grids.r.N, grids.θ.N, grids.ζ.N)
        end
        
    end

    return ϕ, ϕft

end


"""
    allocate_phi_arrays(grids::FFFGridsT, nevals=0; deriv=false)

Allocates empty arrays used for computing phi and its fourier transformation. If nevals is non-zero creates the global phi structures that store all eigenvalues.
"""
function allocate_phi_arrays(grids::FFFGridsT, nevals=0; deriv=false)

    if deriv

        if nevals == 0
            ϕ = Array{ComplexF64}(undef, grids.r.N, grids.θ.N, grids.ζ.N, 8)
            ϕft = Array{ComplexF64}(undef, grids.r.N, grids.θ.N, grids.ζ.N, 8)
        else
            ϕ = Array{ComplexF64}(undef, nevals, grids.r.N, grids.θ.N, grids.ζ.N, 8)
            ϕft = Array{ComplexF64}(undef, nevals, grids.r.N, grids.θ.N, grids.ζ.N, 8)
        end
        
    else
        if nevals == 0
            ϕ = Array{ComplexF64}(undef, grids.r.N, grids.θ.N, grids.ζ.N)
            ϕft = Array{ComplexF64}(undef, grids.r.N, grids.θ.N, grids.ζ.N)
        else
            ϕ = Array{ComplexF64}(undef, nevals, grids.r.N, grids.θ.N, grids.ζ.N)
            ϕft = Array{ComplexF64}(undef, nevals, grids.r.N, grids.θ.N, grids.ζ.N)
        end
        
    end

    return ϕ, ϕft


end


end