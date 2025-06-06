"""
    reconstruct_phi!(efunc::Array{ComplexF64}, grids::FSSGridsT, ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, plan::FFTW.FFTWPlan)

Reconstructs the 1d eigenfunction output back to the 3d grid.
"""
function reconstruct_phi!(efunc::Array{ComplexF64}, grids::FSSGridsT, ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, plan::FFTW.FFTWPlan)


    for i in 1:2:matrix_size(grids)

        #note these are the indicies.
        x1, x2, x3, _ = index_to_grid(i, grids)
        ϕft[x1, x2, x3] = efunc[i]

    end
    ft_phi!(ϕ, ϕft, grids, plan)
end


"""
    reconstruct_phi!(efunc::Array{ComplexF64}, grids::FFSGridsT, ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, plan::FFTW.FFTWPlan)

Reconstructs the 1d eigenfunction output back to the 3d grid.
"""
function reconstruct_phi!(efunc::Array{ComplexF64}, grids::FFSGridsT, ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, plan::FFTW.FFTWPlan)


    m = grids.x2.pf

    x2grid = periodic_grid(grids.x2)

    for i in 1:4:matrix_size(grids)

        #note these are the indicies.
        x1, x2, x3, _ = index_to_grid(i, grids)
        ϕft[x1, x2, x3] = efunc[i] * exp(1im * m * x2grid[x2])

    end

    ft_phi!(ϕ, ϕft, grids, plan)
end



"""
    reconstruct_phi!(efunc::Array{ComplexF64}, grids::FFFGridsT, ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, plan::FFTW.FFTWPlan)

Reconstructs the 1d eigenfunction output back to the 3d grid.
"""
function reconstruct_phi!(efunc::Array{ComplexF64}, grids::FFFGridsT, ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, plan::FFTW.FFTWPlan)


    m = grids.x2.pf
    n = grids.x3.pf

    _, x2grid, x3grid = inst_grids(grids)

    for i in 1:8:matrix_size(grids)

        #note these are the indicies.
        x1, x2, x3, _ = index_to_grid(i, grids)
        ϕ[x1, x2, x3] = efunc[i] * exp(1im * (m * x2grid[x2] + n * x3grid[x3]))

    end
    ft_phi!(ϕ, ϕft, grids, plan)
end


"""
    reconstruct_phi!(efunc::Array{ComplexF64}, grids::FSSGridsT, ϕ::Array{ComplexF64, 4}, ϕft::Array{ComplexF64, 4}, plan::FFTW.FFTWPlan)

Reconstructs the 1d eigenfunction output back to the 3d grid.
"""
function reconstruct_phi!(efunc::Array{ComplexF64}, grids::FSSGridsT, ϕ::Array{ComplexF64, 4}, ϕft::Array{ComplexF64, 4}, plan::FFTW.FFTWPlan)


    for i in 1:1:matrix_size(grids)

        #note these are the indicies.
        x1, x2, x3, hs = index_to_grid(i, grids)
        ϕft[x1, x2, x3, hs] = efunc[i]

    end
    ft_phi!(ϕ, ϕft, grids, plan)
end


"""
    reconstruct_phi!(efunc::Array{ComplexF64}, grids::FFSGridsT, ϕ::Array{ComplexF64, 4}, ϕft::Array{ComplexF64, 4}, plan::FFTW.FFTWPlan)

Reconstructs the 1d eigenfunction output back to the 3d grid.
"""
function reconstruct_phi!(efunc::Array{ComplexF64}, grids::FFSGridsT, ϕ::Array{ComplexF64, 4}, ϕft::Array{ComplexF64, 4}, plan::FFTW.FFTWPlan)


    m = grids.x2.pf

    x2grid = periodic_grid(grids.x2)

    for i in 1:matrix_size(grids)

        #note these are the indicies.
        x1, x2, x3, hs = index_to_grid(i, grids)
        ϕft[x1, x2, x3, hs] = efunc[i] * exp(1im * m * x2grid[x2])

    end

    ft_phi!(ϕ, ϕft, grids, plan)
end



"""
    reconstruct_phi!(efunc::Array{ComplexF64}, grids::FFFGridsT, ϕ::Array{ComplexF64, 4}, ϕft::Array{ComplexF64, 4}, plan::FFTW.FFTWPlan)

Reconstructs the 1d eigenfunction output back to the 3d grid.
"""
function reconstruct_phi!(efunc::Array{ComplexF64}, grids::FFFGridsT, ϕ::Array{ComplexF64, 4}, ϕft::Array{ComplexF64, 4}, plan::FFTW.FFTWPlan)


    m = grids.x2.pf
    n = grids.x3.pf

    _, x2grid, x3grid = inst_grids(grids)

    for i in 1:matrix_size(grids)

        #note these are the indicies.
        x1, x2, x3, hs = index_to_grid(i, grids)
        ϕ[x1, x2, x3, hs] = efunc[i] * exp(1im * (m * x2grid[x2] + n * x3grid[x3]))

    end
    ft_phi!(ϕ, ϕft, grids, plan)
end
