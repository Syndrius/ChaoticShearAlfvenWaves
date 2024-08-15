
"""
    compute_boundary_inds(grids::FFSGridsT)

Computes the indicies of the matrices that correspond to the radial boundaries. Boundaries are ϕ(a) = ϕ(0) = 0. Periodic boundaries are handled elsewhere.
"""
function compute_boundary_inds(grids::FFSGridsT)
    
    Nr = grids.r.N
    Nθ = grids.θ.N
    Nn = grids.ζ.count
    #note this only gets r boundaries, θ is assumed periodic and handles elsewhere.
    
    #unsure about this as ever...
    #may need to do grid_to_ind first...
    left_boundary1 = 1:4:4*Nθ * Nn
    left_boundary2 = 2:4:4*Nθ * Nn
    
    #derivs is also set to zero for lfr
    #not sure about this tbh.
    #left_boundary3 = 3:4:4*Nθ * Nn
   
    right_boundary1 = 1 + (Nr - 1) * 4 * Nθ * Nn:4:4*Nr * Nθ * Nn
    right_boundary2 = 2 + (Nr - 1) * 4 * Nθ * Nn:4:4*Nr * Nθ * Nn

    return vcat(left_boundary1, left_boundary2, right_boundary1, right_boundary2)
    #return vcat(left_boundary1, left_boundary2, left_boundary3, right_boundary1, right_boundary2)

end



"""
    grid_to_index(rind::Int64, θind::Int64, ζind::Int64, hr::Int64, hθ::Int64, grids::FFSGridsT)

Converts the index of our 3d grid to the appropriate place in the matrix.
"""
function grid_to_index(rind::Int64, θind::Int64, ζind::Int64, hr::Int64, hθ::Int64, grids::FFSGridsT)

    Nθ = grids.θ.N
    Nn = grids.ζ.count
    
    #matrix is structured by segments following the grid.
    #Follows the pattern (0, 0, 0) -> (0, 0, 1) -> ... -> (0, 0, Nn) -> (0, 1, 0) -> (0, Nθ, Nn) -> (1, 0, 0) etc.

    #Each segment contains the 4 basis function combinations at each node. So each segment is made form 4 indicies.
    #rh1θh1ζ, rh1θ2ζ, rh2θh1ζ, rh2θh2ζ.

    #(subtract 1 due to julia indexing starting from 1)
    # mod is for periodicity, so that if θind==Nθ, it is mapped to θind = 1.

    #TODO describe this

    #first we find the appropriate segment for our indicies.
    # - Each increment of ζind moves us one segment.
    # - Each increment of θind moves us Nζ segments.
    # - Each increment of rind moves us Nθ * Nζ segments.
    segment = (ζind-1) + Nn * mod(θind-1 + grid_id[hθ], Nθ) + Nθ * Nn * (rind - 1 + grid_id[hr]) 

    #Once we have the segment, we need to find the correct index.
    return 1 + basis_id[hθ] + 2 * basis_id[hr] + 4 * segment


end



"""
    reconstruct_phi(efuncs::Array{ComplexF64, 2}, nevals::Int64, grids::FFSGridsT)

Reconstructs the 1d eigenfunction output back to the 3d grid.
"""
function reconstruct_phi(efuncs::Array{ComplexF64, 2}, nevals::Int64, grids::FFSGridsT)
    phi = zeros(ComplexF64, nevals, grids.r.N, grids.θ.N, grids.ζ.count)
    #maybe one day we will want dphidr???
    #note that this is the same for both!

    #m = grids.θ.pf

    #_, θgrid, _, _, _ = instantiate_grids(grids)

    for i in 1:4:matrix_dim(grids)

        #note these are the indicies.
        r, θ, ζ, hs = index_to_grid(i, grids)
        phi[:, r, θ, ζ] = efuncs[i, :]

        #if hs == 1
            #may be the wrong way around!
            #this doesn't seem to have worked as expected tbh!
            #hard to tell what is needed here!
            #phi[:, r, θ, ζ] = efuncs[i, :] .* exp(1im * m * θgrid[θ])
        #end
    end
    return phi
end


"""
    index_to_grid(i::Int64, grids::FFSGridsT)

Reverses the transformation from 3d grid to matrix, i.e. maps matrix index to appropriate point in the grid, used for plotting.
"""
function index_to_grid(i::Int64, grids::FFSGridsT)
    r = div(i-1, 4*grids.θ.N*grids.ζ.count) + 1
    θ = mod(div(i-1, 4*grids.ζ.count), grids.θ.N) + 1
    ζ = mod(div(i-1, 4), grids.ζ.count) + 1

    h = mod(i-1, 4) + 1 #ie even or odd

    return r, θ, ζ, h

end




"""
    function local_to_global(rnode::Int64, θnode::Int64, ξr::Array{Float64}, ξθ::Array{Float64}, rgrid::Array{Float64}, θgrid::StepRangeLen)

Converts the local grid where the finite elements are defined to the global coordinates.
"""
function local_to_global(rnode::Int64, θnode::Int64, ξr::Array{Float64}, ξθ::Array{Float64}, rgrid::Array{Float64}, θgrid::StepRangeLen)

    #hopefully handles periodicity!
    #this probably assumes θgrid is evenly spaced, at least across the periodic part.
    θnode = mod(θnode, length(θgrid)-1) + 1

    dr = rgrid[rnode+1] - rgrid[rnode]
    dθ = θgrid[θnode+1] - θgrid[θnode]

    #first we map the (-1, 1) local coords to (0, 1)
    mpr = @. (ξr+1) / 2
    mpθ = @. (ξθ+1) / 2

    #then we sale by the width of grid points
    mpr = mpr .* dr
    mpθ = mpθ .* dθ

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    rglobal = mpr .+ rgrid[rnode]
    θglobal = mpθ .+ θgrid[θnode]

    
    return rglobal, θglobal, dr, dθ
end 



"""
    matrix_dim(grids::FFSGridsT)

Computes the dimension of the matrix.
"""
function matrix_dim(grids::FFSGridsT)

    return 4 * grids.r.N * grids.θ.N * grids.ζ.count
end



###############################################################
#Legacy function! for when basis functions where combined.

#we can probbaly still do this with segment formalism, but it may not be as simple.
#unfort has the same args as other case so will need a different name.
function grid_to_index_comb(rind::Int64, θind::Int64, ζind::Int64, hind::Int64, grids::FFSGridsT)

    Nθ = grids.θ.N
    Nn = grids.ζ.count
    #id's for the heremite basis
    #tells us which node the basis belongs to, 0 for left, 1 for right of current element
    #these should probably not be redefined all the time...
    
    #tells us which of the basis vectors is being considered here
    #0 reps the first which has restrictions on the value at the boundaries
    #1 reps the second which has restrictions on the derivative at the boundaries
    

    #matrix will be structured so [r1θ1ζ1, r1θ1ζ2, ...r1θ1ζnn, r1θ2ζ1, ..., r1θnmζnn, r2θ1ζ1, .., rnr, θnm, ζnn]

    #but each point has to have four values, for each 2d hermite basis.
    #where r_11 is first radial point, first Hermite, r_12 is first radial second Hermite ie the derivative one
    #so it will actually be [r11θ1ζ1, r12θ1ζ1, r11θ1mζ2...]

    #so it will actually be [r11θ11ζ1, r12θ11ζ1, r11θ12ζ1, r12θ12ζ1, r11θ1ζ2...]

    #matrix goes [r1θ1ζ1, r1θ1ζ2, ..., r1θ1ζn, r1θ2ζ1, ..., r1θnζn, r2θ1ζ1, ...]

    #with each segment being, #r1θ1ζ1, made up of four case
    #[rh1θh1ζ, rh1θh2ζ, rh2θh1ζ, rh2θh2ζ]

    #start at 1 for julia arrays
    #shift across 1 if we are dealing with the derivative basis of r
    #shift across 2 if we are dealing with the derivative basis of θ
    #i.e. basis funcs goes S(0, 0) = 1, dSr(0, 0) = 1, dSθ(0, 0) = 1, ddSrθ(0, 0) = 1, where these derivs refer to basis functions maxima!
    #shift by 4 for each ζind, ie the number of basis funcs at each node.
    #shift by 4 * ζpoints for each θ, i.e the number of basis funcs for each ζ node.
    #shift by 4 * θpoints * ζpoints for each r, i.e number of basis funcs for each θ and ζ node.
    #so we run across ζ values, then increment θ, repeat until finish θ, then increment r.
    #next we shift the θind up one if our S is on the neighbouring node
    #same for r.
    #mod θ maps last θ point back to same index as first point 

    #need to do periodicty!
    return 1 + basis_id_θ[hind] + 2 * basis_id_r[hind] + 4 * (ζind-1) + 4 * Nn * mod(θind-1 + grid_id_θ[hind], Nθ) + 4 * Nθ * Nn * (rind - 1 + grid_id_r[hind]) 



    #the index from 1 is going to change things a bit.
    #can be made clearer but does seem to be working ish.
    #return 1+4*Nm*Nn*(rind-1) + basis_id[sind] + 2*nn*(θind-1) + 2*(ζind-1) + 2*nm*nn*grid_id[hbind]

end