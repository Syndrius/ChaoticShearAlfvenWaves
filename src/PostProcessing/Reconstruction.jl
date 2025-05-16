"""
    reconstruct_phi(efuncs::Array{ComplexF64, 2}, nevals::Int64, grids::FSSGridsT)

Reconstructs the 1d eigenfunction output back to the 3d grid. Single efunc case.
"""
function reconstruct_phi!(efunc::Array{ComplexF64}, grids::FSSGridsT, ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, plan::FFTW.FFTWPlan)


    for i in 1:2:matrix_size(grids)

        #note these are the indicies.
        x1, x2, x3, hs = index_to_grid(i, grids)
        ϕft[x1, x2, x3] = efunc[i]

    end

    ft_phi!(ϕ, ϕft, grids, plan)
end



"""
    reconstruct_phi(efuncs::Array{ComplexF64, 2}, nevals::Int64, grids::FFSGridsT)

Reconstructs the 1d eigenfunction output back to the 3d grid. This case for a single eigenfunction.
"""
function reconstruct_phi!(efunc::Array{ComplexF64}, grids::FFSGridsT, ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, plan::FFTW.FFTWPlan)


    m = grids.x2.pf

    x2grid = periodic_grid(grids.x2)

    for i in 1:4:matrix_size(grids)

        #note these are the indicies.
        x1, x2, x3, hs = index_to_grid(i, grids)
        ϕft[x1, x2, x3] = efunc[i] * exp(1im * m * x2grid[x2])

    end
    ft_phi!(ϕ, ϕft, grids, plan)
end




"""
    reconstruct_phi(efuncs::Array{ComplexF64}, grids::FFFGridsT)

Reconstructs the 1d eigenfunction output back to the 3d grid. Only a single efunc.
"""
function reconstruct_phi!(efunc::Array{ComplexF64}, grids::FFFGridsT, ϕ::Array{ComplexF64, 3}, ϕft::Array{ComplexF64, 3}, plan::FFTW.FFTWPlan)


    m = grids.x2.pf
    n = grids.x3.pf

    _, x2grid, x3grid = inst_grids(grids)

    for i in 1:8:matrix_size(grids)

        #note these are the indicies.
        x1, x2, x3, hs = index_to_grid(i, grids)
        ϕ[x1, x2, x3] = efunc[i] * exp(1im * (m * x2grid[x2] + n * x3grid[x3]))

    end
    ft_phi!(ϕ, ϕft, grids, plan)
end



#######################################################
#TODO other cases including the derivatives 

"""
    reconstruct_phi(efuncs::Array{ComplexF64, 2}, nevals::Int64, grids::FFFGridsT)

Reconstructs the 1d eigenfunction output back to the 3d grid.
"""
function reconstruct_phi!(efunc::Array{ComplexF64}, grids::FFFGridsT, ϕ::Array{ComplexF64, 4})
    #phi = zeros(ComplexF64, grids.r.N, grids.x2.N, grids.x3.N, 8)
    #phi = Array{ComplexF64}(undef, grids.r.N, grids.x2.N, grids.x3.N, 8)
    #maybe one day we will want dphidr???
    #note that this is the same for both!

    #TODO

    m = grids.x2.pf
    n = grids.x3.pf

    _, x2grid, x3grid = inst_grids(grids)

    #should skip over derivative inds.
    for i in 1:1:matrix_dim(grids)

        

        #note these are the indicies.
        r, x2, x3, hs = index_to_grid(i, grids)

        
        #the 8 should fix this
        ϕ[r, x2, x3, hs] = efunc[i] * exp(1im * (m * x2grid[x2] + n * x3grid[x3]))

        #
            #may be the wrong way around!
            #this doesn't seem to have worked as expected tbh!
            #phi[:, r, x2, x3] = efuncs[i, :] #.* exp(1im * (m * x2grid[x2] + n * x3grid[x3]))
        #end
    end

end





"""
    reconstruct_phi(efuncs::Array{ComplexF64, 2}, nevals::Int64, grids::FFSGridsT)

Reconstructs the 1d eigenfunction output back to the 3d grid. This case for a single eigenfunction.
"""
function reconstruct_phi!(efunc::Array{ComplexF64}, grids::FFSGridsT)

    phi = zeros(ComplexF64, grids.r.N, grids.x2.N, grids.x3.count, 4)
    #maybe one day we will want dphidr???
    #note that this is the same for both!

    m = grids.x2.pf

    _, x2grid, _, _, _ = instantiate_grids(grids)

    for i in 1:1:matrix_dim(grids)

        #note these are the indicies.
        x1, x2, x3, hs = index_to_grid(i, grids)
        phi[x1, x2, x3, hs] = efunc[i] * exp(1im * m * x2grid[x2])

        #if hs == 1
            #may be the wrong way around!
            #this doesn't seem to have worked as expected tbh!
            #hard to tell what is needed here!
            #phi[:, r, x2, x3] = efuncs[i, :] .* exp(1im * m * x2grid[x2])
        #end
    end
    return phi
end



