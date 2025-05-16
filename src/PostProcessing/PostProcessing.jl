"""

This modules handles the processing of the eigenvalues and eigenfunctions returned by the solver into convenient forms. In particular, the eigenvalues are normalised to the eigenfrequency at the axis, and the eigenfunctions are converted from 1d arrays into 3d arrays reflecting the grid. Additionally, the eigenfunctions are also Fourier transformed in x2 and x3 so the mode structure can be viewed.
"""
module PostProcessing

using FFTW
using EllipsisNotation

using ..Structures
using ..Indexing
using ..Basis


export post_process


include("Reconstruction.jl")


include("FT.jl")


include("Continuum.jl")



"""
    post_process(evals::AbstractArray, efuncs::Array{ComplexF64}, grids::GridsT, geo::GeoParamsT, deriv=false::Bool)

Processes the eigenvalues and eigenfunctions into useful forms. This involves recosntructing the 1d eigenfunction into the 3d grid representation for ϕ, and taking the fourier transform in x2 and x3. The fourier transform is used to label each eigenvalue with its dominant mode and radial location. The eigenvalues are also normalised to the Alfven frequency.
"""
function post_process(evals::AbstractArray, efuncs::Array{ComplexF64}, grids::GridsT, geo::GeoParamsT, deriv=false::Bool)

    #deriv is only really valid for fff for now, and probably for a while lol.

    #allocates global arrays for storing everyeigenvalue
    ϕ_g, ϕft_g = allocate_phi_arrays(grids, length(evals), deriv=deriv)
    #allocates placeholder arrays used during computation.
    ϕp, ϕpft = allocate_phi_arrays(grids, deriv=deriv)

    #array to store the radial value of each eigenvalue.
    x1ms = Array{Float64}(undef, length(evals))

    #array to store the mode labels
    mode_labs = Tuple{Int, Int}[] 

    #create a plan for fourier transforming.
    plan = create_ft_plan(ϕpft, grids)

    x1grid = inst_grids(grids)[1]

    #arrays to store the maximum value of ϕft and the correspondning r value.
    x1marray = Array{Int64}(undef, grids.x2.N, grids.x3.N)
    ϕmarray = Array{Float64}(undef, grids.x2.N, grids.x3.N)

    
    for i in 1:length(evals)

        #converts the eigenfunction back into the 3d ϕ, including its fourier transformation in x2 and x3
        reconstruct_phi!(efuncs[:, i], grids, ϕp, ϕpft, plan)

        #finds the dominant mode and its radial location for continuum reconstruction.
        x1ind, mode_lab = label_mode(ϕpft, grids, x1marray, ϕmarray)

        ϕ_g[i, ..] = ϕp
        ϕft_g[i, ..] = ϕpft

        x1ms[i] = x1grid[x1ind]
        push!(mode_labs, mode_lab)
    end

    #should only take abs for ideal cases, this is just to avoid
    #tiny negatives that should be zero
    #alterantively, as per cka, just set any tiny negatives to zero.
    ω = geo.R0 .* sqrt.(abs.(evals))
    #ω = geo.R0 .* sqrt.(evals)

    evals = EvalsT(ω, x1ms, mode_labs)

    return evals, ϕ_g, ϕft_g

end


"""
    allocate_phi_arrays(grids::FSSGridsT, nevals=0; deriv=false)

Allocates empty arrays used for computing phi and its fourier transformation. If nevals is non-zero creates the global phi structures that store all eigenvalues.
"""
function allocate_phi_arrays(grids::FSSGridsT, nevals=0; deriv=false)

    #extends the domain for ifft cases for smoother plotting.
    x2size = ifft_size(grids.x2)
    x3size = ifft_size(grids.x3)

    if deriv

        if nevals == 0
            ϕ = Array{ComplexF64}(undef, grids.x1.N, x2size, x3size, 8)
            ϕft = Array{ComplexF64}(undef, grids.x1.N, grids.x2.N, grids.x3.N, 8)
        else
            ϕ = Array{ComplexF64}(undef, nevals, grids.x1.N, x2size, x3size, 8)
            ϕft = Array{ComplexF64}(undef, nevals, grids.x1.N, grids.x2.N, grids.x3.N, 8)
        end
        #placeholder memory for each eigenfunction.
        
    else
        if nevals == 0
            ϕ = Array{ComplexF64}(undef, grids.x1.N, x2size, x3size)
            ϕft = Array{ComplexF64}(undef, grids.x1.N, grids.x2.N, grids.x3.N)
        else
            ϕ = Array{ComplexF64}(undef, nevals, grids.x1.N, x2size, x3size)
            ϕft = Array{ComplexF64}(undef, nevals, grids.x1.N, grids.x2.N, grids.x3.N)
        end
        
    end

    return ϕ, ϕft

end



"""
    allocate_phi_arrays(grids::FFSGridsT, nevals=0; deriv=false)

Allocates empty arrays used for computing phi and its fourier transformation. If nevals is non-zero creates the global phi structures that store all eigenvalues.
"""
function allocate_phi_arrays(grids::FFSGridsT, nevals=0; deriv=false)

    x3size = ifft_size(grids.x3)

    if deriv

        if nevals == 0
            ϕ = Array{ComplexF64}(undef, grids.x1.N, grids.x2.N, x3size, 8)
            ϕft = Array{ComplexF64}(undef, grids.x1.N, grids.x2.N, grids.x3.N, 8)
        else
            ϕ = Array{ComplexF64}(undef, nevals, grids.x1.N, grids.x2.N, x3size, 8)
            ϕft = Array{ComplexF64}(undef, nevals, grids.x1.N, grids.x2.N, grids.x3.N, 8)
        end
        
    else
        if nevals == 0
            ϕ = Array{ComplexF64}(undef, grids.x1.N, grids.x2.N, x3size)
            ϕft = Array{ComplexF64}(undef, grids.x1.N, grids.x2.N, grids.x3.N)
        else
            ϕ = Array{ComplexF64}(undef, nevals, grids.x1.N, grids.x2.N, x3size)
            ϕft = Array{ComplexF64}(undef, nevals, grids.x1.N, grids.x2.N, grids.x3.N)
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
            ϕ = Array{ComplexF64}(undef, grids.x1.N, grids.x2.N, grids.x3.N, 8)
            ϕft = Array{ComplexF64}(undef, grids.x1.N, grids.x2.N, grids.x3.N, 8)
        else
            ϕ = Array{ComplexF64}(undef, nevals, grids.x1.N, grids.x2.N, grids.x3.N, 8)
            ϕft = Array{ComplexF64}(undef, nevals, grids.x1.N, grids.x2.N, grids.x3.N, 8)
        end
        
    else
        if nevals == 0
            ϕ = Array{ComplexF64}(undef, grids.x1.N, grids.x2.N, grids.x3.N)
            ϕft = Array{ComplexF64}(undef, grids.x1.N, grids.x2.N, grids.x3.N)
        else
            ϕ = Array{ComplexF64}(undef, nevals, grids.x1.N, grids.x2.N, grids.x3.N)
            ϕft = Array{ComplexF64}(undef, nevals, grids.x1.N, grids.x2.N, grids.x3.N)
        end
        
    end

    return ϕ, ϕft


end


end
