
"""
    compute_boundary_inds(grids::FFFGridsT)

Computes the indicies of the matrices that correspond to the radial boundaries. Boundaries are ϕ(a) = ϕ(0) = 0. Periodic boundaries are handled elsewhere.
"""
function compute_boundary_inds(grids::FFFGridsT)
    
    Nr = grids.r.N
    Nθ = grids.θ.N
    Nζ = grids.ζ.N
    #note this only gets r boundaries, θ is assumed periodic and handles elsewhere.
    
    #unsure about this as ever...
    #may need to do grid_to_ind first...
    left_boundary1 = 1:8:8*Nθ * Nζ
    left_boundary2 = 2:8:8*Nθ * Nζ
    left_boundary3 = 3:8:8*Nθ * Nζ
    left_boundary4 = 4:8:8*Nθ * Nζ
   
    right_boundary1 = 1 + (Nr - 1) * 8 * Nθ * Nζ:8:8*Nr * Nθ * Nζ
    right_boundary2 = 2 + (Nr - 1) * 8 * Nθ * Nζ:8:8*Nr * Nθ * Nζ
    right_boundary3 = 3 + (Nr - 1) * 8 * Nθ * Nζ:8:8*Nr * Nθ * Nζ
    right_boundary4 = 4 + (Nr - 1) * 8 * Nθ * Nζ:8:8*Nr * Nθ * Nζ

    return vcat(left_boundary1, left_boundary2, left_boundary3, left_boundary4, right_boundary1, right_boundary2, right_boundary3, right_boundary4)

end


"""
    grid_to_index(rind::Int64, θind::Int64, ζind::Int64, hr::Int64, hθ::Int64, hζ::Int64, grids::FFFGridsT)

Converts the index of our 3d grid to the appropriate place in the matrix.
"""
function grid_to_index(rind::Int64, θind::Int64, ζind::Int64, hr::Int64, hθ::Int64, hζ::Int64, grids::FFFGridsT)

    Nθ = grids.θ.N
    Nζ = grids.ζ.N
    
    #matrix is structured by segments following the grid.
    #Follows the pattern (0, 0, 0) -> (0, 0, 1) -> ... -> (0, 0, Nn) -> (0, 1, 0) -> (0, Nθ, Nn) -> (1, 0, 0) etc.

    #Each segment contains the 8 basis function combinations at each node. So each segment is made form 8 indicies.
    #rh1θh1ζh1, rh1θh1ζh2, rh1θ2ζh1, rh1θh2ζh2, rh2θh1ζh1, rh2θh1ζh2, rh2θh2ζh1, rh2θh2ζh2.

    #(subtract 1 due to julia indexing starting from 1)
    # mod is for periodicity, so that if ζind==Nζ, it is mapped to ζind = 1.


    #first we find the appropriate segment for our indicies.
    # - Each increment of ζind moves us one segment.
    # - Each increment of θind moves us Nζ segments.
    # - Each increment of rind moves us Nθ * Nζ segments.
    segment = mod(ζind-1 + grid_id[hζ], Nζ) + Nζ * mod(θind-1 + grid_id[hθ], Nθ) + Nθ * Nζ * (rind - 1 + grid_id[hr]) 

    #Once we have the segment, we need to find the correct index.
    return 1 + 8 * segment + basis_id[hζ] + 2 * basis_id[hθ] + 4 * basis_id[hr]
end



"""
    local_to_global(rnode::Int64, θnode::Int64, ζnode::Int64, ξr::Array{Float64}, ξθ::Array{Float64}, ξζ::Array{Float64}, rgrid::Array{Float64}, θgrid::StepRangeLen, ζgrid::StepRangeLen)

Converts the local grid where the finite elements are defined to the global coordinates.
"""
function local_to_global(rnode::Int64, θnode::Int64, ζnode::Int64, ξr::Array{Float64}, ξθ::Array{Float64}, ξζ::Array{Float64}, rgrid::Array{Float64}, θgrid::StepRangeLen, ζgrid::StepRangeLen)

    #hopefully handles periodicity!
    #this probably assumes θgrid is evenly spaced, at least across the periodic part.
    #θnode = mod(θnode, length(θgrid)-1) + 1
    #ζnode = mod(ζnode, length(ζgrid)-1) + 1

    if θnode == length(θgrid)
        dθ = 2π + θgrid[1] - θgrid[θnode]
    else
        dθ = θgrid[θnode+1] - θgrid[θnode]
    end

    if ζnode == length(ζgrid)
        dζ = 2π + ζgrid[1] - ζgrid[ζnode]
    else
        dζ = ζgrid[ζnode+1] - ζgrid[ζnode]
    end

    dr = rgrid[rnode+1] - rgrid[rnode]
    

    #first we map the (-1, 1) local coords to (0, 1)
    mpr = @. (ξr+1) / 2
    mpθ = @. (ξθ+1) / 2
    mpζ = @. (ξζ+1) / 2

    #then we sale by the width of grid points
    mpr = mpr .* dr
    mpθ = mpθ .* dθ
    mpζ = mpζ .* dζ

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    rglobal = mpr .+ rgrid[rnode]
    θglobal = mpθ .+ θgrid[θnode]
    ζglobal = mpζ .+ ζgrid[ζnode]

    
    return rglobal, θglobal, ζglobal, dr, dθ, dζ
end 


"""
    index_to_grid(i::Int64, grids::FFFGridsT)

Reverses the transformation from 3d grid to matrix, i.e. maps matrix index to appropriate point in the grid, used for plotting.
"""
function index_to_grid(i::Int64, grids::FFFGridsT)
    r = div(i-1, 8*grids.θ.N*grids.ζ.N) + 1
    θ = mod(div(i-1, 8*grids.ζ.N), grids.θ.N) + 1
    ζ = mod(div(i-1, 8), grids.ζ.N) + 1

    h = mod(i-1, 8) + 1 #ie even or odd

    return r, θ, ζ, h

end


"""
    matrix_dim(grids::FFFGridsT)

Computes the dimension of the matrix.
"""
function matrix_dim(grids::FFFGridsT)

    return 8 * grids.r.N * grids.θ.N * grids.ζ.N
end


"""
    reconstruct_phi(efuncs::Array{ComplexF64, 2}, nevals::Int64, grids::FFFGridsT)

Reconstructs the 1d eigenfunction output back to the 3d grid.
"""
function reconstruct_phi(efuncs::Array{ComplexF64, 2}, nevals::Int64, grids::FFFGridsT)
    phi = zeros(ComplexF64, nevals, grids.r.N, grids.θ.N, grids.ζ.N)
    #maybe one day we will want dphidr???
    #note that this is the same for both!

    #m = grids.θ.pf
    #n = grids.ζ.pf

    #_, θgrid, ζgrid = instantiate_grids(grids)

    #should skip over derivative inds.
    for i in 1:8:matrix_dim(grids)

        #note these are the indicies.
        r, θ, ζ, hs = index_to_grid(i, grids)
        #the 8 should fix this
        phi[:, r, θ, ζ] = efuncs[i, :]

        #if hs == 1
            #may be the wrong way around!
            #this doesn't seem to have worked as expected tbh!
            #phi[:, r, θ, ζ] = efuncs[i, :] #.* exp(1im * (m * θgrid[θ] + n * ζgrid[ζ]))
        #end
    end
    return phi
end