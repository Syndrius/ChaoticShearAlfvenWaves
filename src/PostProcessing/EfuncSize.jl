"""
    allocate_phi_arrays(grids::FSSGridsT, nevals::Int64=0; deriv::Bool=false)

Allocates empty arrays used for computing phi and its fourier transformation. If nevals is non-zero creates the global phi structures that store all eigenvalues.
"""
function allocate_phi_arrays(grids::FSSGridsT, nevals::Int64=0; deriv::Bool=false)

    #extends the domain for ifft cases for smoother plotting.
    x2size = ifft_size(grids.x2)
    x3size = ifft_size(grids.x3)

    if deriv

        if nevals == 0
            #these being 8 instead of 2 is fine as hermite only interacts with the first 2,
            #and this will allow us to maybe add fourier derivatives later.
            ϕ = Array{ComplexF64}(undef, grids.x1.N, x2size, x3size, 8)
            ϕft = Array{ComplexF64}(undef, grids.x1.N, grids.x2.N, grids.x3.N, 8)
        else
            ϕ = Array{ComplexF64}(undef, nevals, grids.x1.N, x2size, x3size, 8)
            ϕft = Array{ComplexF64}(undef, nevals, grids.x1.N, grids.x2.N, grids.x3.N, 8)
        end
        
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
    allocate_phi_arrays(grids::FFSGridsT, nevals::Int64=0; deriv::Bool=false)

Allocates empty arrays used for computing phi and its fourier transformation. If nevals is non-zero creates the global phi structures that store all eigenvalues.
"""
function allocate_phi_arrays(grids::FFSGridsT, nevals::Int64=0; deriv::Bool=false)

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
    allocate_phi_arrays(grids::FFFGridsT, nevals::Int64=0; deriv::Bool=false)

Allocates empty arrays used for computing phi and its fourier transformation. If nevals is non-zero creates the global phi structures that store all eigenvalues.
"""
function allocate_phi_arrays(grids::FFFGridsT, nevals::Int64=0; deriv::Bool=false)

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

