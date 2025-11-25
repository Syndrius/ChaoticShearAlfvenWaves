"""
    create_ft_plan(ϕ::Array{ComplexF64}, grids::FSSGridsT)

Creates a fourier transform plan to speed up repeated ft. For FSS this is not actually used, but keeps type information.
"""
function create_ft_plan(ϕ::Array{ComplexF64}, grids::FSSGridsT)
    return plan_fft(ϕ) #not actually used, needed for type information.
end


"""
    create_ft_plan(ϕ::Array{ComplexF64}, grids::FFSGridsT)

Creates a fourier transform plan to speed up repeated ft. 
"""
function create_ft_plan(ϕ::Array{ComplexF64}, grids::FFSGridsT)
    return plan_fft(ϕ, [2])
end


"""
    create_ft_plan(ϕ::Array{ComplexF64}, grids::FFFGridsT)

Creates a fourier transform plan to speed up repeated ft. 
"""
function create_ft_plan(ϕ::Array{ComplexF64}, grids::FFFGridsT)
    return plan_fft(ϕ, [2, 3])
end


"""
    ft_phi!(ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, grids::FSSGridsT, plan::FFTW.FFTWPlan)

Computes the fourier and non-fourier representation of phi once it has been reconstructed.
"""
function ft_phi!(ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, grids::FSSGridsT, plan::FFTW.FFTWPlan)

    x2grid_size = ifft_size(grids.x2)
    x3grid_size = ifft_size(grids.x3)

    x2grid = periodic_grid(x2grid_size)
    x3grid = periodic_grid(x3grid_size)

    mlist = mode_list(grids.x2)
    nlist = mode_list(grids.x3)

    #here an 'inverse fourier transform' is applied, custom loop as only specific modes are used.
    #note that this can be slow.
    ϕ .= 0
    for i in 1:grids.x1.N
        for j in 1:1:length(mlist)
            for k in 1:1:length(nlist)

                ϕ[i, :, :] += ϕft[i, j, k] .* exp.(1im * mlist[j] .* x2grid .+ 1im * nlist[k] .* x3grid' )
            end
        end
    end

end



"""
    ft_phi!(ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, grids::FSSGridsT, plan::FFTW.FFTWPlan)

Computes the fourier and non-fourier representation of phi once it has been reconstructed.
"""
function ft_phi!(ϕ::Array{ComplexF64, 4}, ϕft::Array{ComplexF64, 4}, grids::FSSGridsT, plan::FFTW.FFTWPlan)

    x2grid_size = ifft_size(grids.x2)
    x3grid_size = ifft_size(grids.x3)

    x2grid = periodic_grid(x2grid_size)
    x3grid = periodic_grid(x3grid_size)

    mlist = mode_list(grids.x2)
    nlist = mode_list(grids.x3)

    #here an 'inverse fourier transform' is applied, custom loop as only specific modes are used.
    ϕ .= 0
    for i in 1:grids.x1.N
        for j in 1:1:length(mlist)
            for k in 1:1:length(nlist)
                for l in 1:2

                    ϕ[i, :, :, l] += ϕft[i, j, k, l] .* exp.(1im * mlist[j] .* x2grid .+ 1im * nlist[k] .* x3grid' )
                end
            end
        end
    end

end



"""
    ft_phi!(ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, grids::FFSGridsT, plan::FFTW.FFTWPlan)

Computes the fourier and non-fourier representation of phi once it has been reconstructed.
"""
function ft_phi!(ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, grids::FFSGridsT, plan::FFTW.FFTWPlan)

    x3grid_size = ifft_size(grids.x3)

    x3grid = periodic_grid(x3grid_size)

    nlist = mode_list(grids.x3)

    #here an 'inverse fourier transform' is applied, custom loop as only specific modes are used.
    ϕ .= 0
    for i in 1:grids.x1.N
        
        for j in 1:grids.x2.N
            for k in 1:1:length(nlist)

                ϕ[i, j, :] += ϕft[i, j, k] .* exp.(1im * nlist[k] .* x3grid)
            end
        end
    end

    #now fourier transform in x2
    ϕft .= plan * ϕft

end



"""
    ft_phi!(ϕ::Array{ComplexF64, 4}, ϕft::Array{ComplexF64, 4}, grids::FFSGridsT, plan::FFTW.FFTWPlan)

Computes the fourier and non-fourier representation of phi once it has been reconstructed.
"""
function ft_phi!(ϕ::Array{ComplexF64, 4}, ϕft::Array{ComplexF64, 4}, grids::FFSGridsT, plan::FFTW.FFTWPlan)

    x3grid_size = ifft_size(grids.x3)

    x3grid = periodic_grid(x3grid_size)

    nlist = mode_list(grids.x3)

    #here an 'inverse fourier transform' is applied, custom loop as only specific modes are used.
    ϕ .= 0
    for i in 1:grids.x1.N
        
        for j in 1:grids.x2.N
            for k in 1:1:length(nlist)

                ϕ[i, j, :, :] += ϕft[i, j, k, :]' .* exp.(1im * nlist[k] .* x3grid)
            end
        end
    end

    #now fourier transform in x2
    ϕft .= plan * ϕft

end


"""
    ft_phi!(ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, grids::FFSGridsT, plan::FFTW.FFTWPlan)

Computes the fourier and non-fourier representation of phi once it has been reconstructed.
"""
function ft_phi!(ϕ::Array{ComplexF64}, ϕft::Array{ComplexF64}, grids::FFFGridsT, plan::FFTW.FFTWPlan)
    
    #FFF case is simply a direct ft of x2 and x3
    #note that this case is the same when the derivative is included
    ϕft .= plan * ϕ

end


