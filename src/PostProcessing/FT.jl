

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

    θgrid_size = ifft_size(grids.θ)
    ζgrid_size = ifft_size(grids.ζ)

    θgrid = periodic_grid(θgrid_size)
    ζgrid = periodic_grid(ζgrid_size)

    mlist = mode_list(grids.θ)
    nlist = mode_list(grids.ζ)

    #here an 'inverse fourier transform' is applied, custom loop as only specific modes are used.
    ϕ .= 0
    for i in 1:grids.r.N
        for j in 1:1:length(mlist)
            for k in 1:1:length(nlist)

                ϕ[i, :, :] += ϕft[i, j, k] .* exp.(1im * mlist[j] .* θgrid .+ 1im * nlist[k] .* ζgrid' )
            end
        end
    end

end



"""
    ft_phi!(ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, grids::FFSGridsT, plan::FFTW.FFTWPlan)

Computes the fourier and non-fourier representation of phi once it has been reconstructed.
"""
function ft_phi!(ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, grids::FFSGridsT, plan::FFTW.FFTWPlan)

    ζgrid_size = ifft_size(grids.ζ)

    ζgrid = periodic_grid(ζgrid_size)

    nlist = mode_list(grids.ζ)

    #here an 'inverse fourier transform' is applied, custom loop as only specific modes are used.
    ϕ .= 0
    for i in 1:grids.r.N
        
        for j in 1:grids.θ.N
            for k in 1:1:length(nlist)

                ϕ[i, j, :] += ϕft[i, j, k] .* exp.(1im * nlist[k] .* ζgrid)
            end
        end
    end

    #now fourier transform in θ
    ϕft .= plan * ϕft

end


"""
    ft_phi!(ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, grids::FFSGridsT, plan::FFTW.FFTWPlan)

Computes the fourier and non-fourier representation of phi once it has been reconstructed.
"""
function ft_phi!(ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, grids::FFFGridsT, plan)
    
    #FFF case is simply a direct ft of θ and ζ
    ϕft .= plan * ϕ

end


#TODO derivtive versions of ft_phi.

