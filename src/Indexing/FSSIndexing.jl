

"""
    compute_boundary_inds(grids::FSSGridsT)

Computes the indicies of the matrices that correspond to the radial boundaries. Boundaries are ϕ(a) = ϕ(0) = 0. Periodic boundaries are handled elsewhere.
"""
function compute_boundary_inds(grids::FSSGridsT)

    Nr = grids.r.N
    Nm = grids.θ.count
    Nn = grids.ζ.count

    #_, mlist, _ = Inputs.sm_grid(grids.θ) #think this is pointless!

    #for radiative case the deriv needs to be zero as well.
    left_boundary = 1:2:2*Nm*Nn

    #flr stuff.
    #left_boundary = 1:1:2*Nm*Nn
    #tried to make the distinction between nθ and nm clearer.
    #think this is a waste of time, hard to know for sure though!
    """
    if 0 in mlist
        #fu and berk radiative paper seem to imply it should still be m=1 even for phi case...
        zero_ind = argmin(abs.(mlist))
        left_boundary = collect(left_boundary)
        #more complicated regularization condition
        #shift the m=0 conditions to the derivative,
        for i in (length(mlist) - zero_ind)*Nn+1:(zero_ind)*Nn
            left_boundary[i] += 1
        end
    end
    """

    #could probably use grid_to_index for this
    right_boundary = 1+(Nr-1)*2*Nm*Nn:2:Nr*2*Nm*Nn

    return vcat(left_boundary, right_boundary)
end



"""
    grid_to_index(rind::Int64, θind::Int64, ζind::Int64, hbind::Int64, grids::FSSGridsT)

Converts the index of our 3d grid to the appropriate place in the matrix.
"""
function grid_to_index(rind::Int64, θind::Int64, ζind::Int64, hr::Int64, grids::FSSGridsT)

    Nm = grids.θ.count
    Nn = grids.ζ.count

    #id's for the heremite basis
    #tells us which node the basis belongs to, 0 for left, 1 for right of current element
    #these should probably not be redefined all the time...
    
    #tells us which of the basis vectors is being considered here
    #0 reps the first which has restrictions on the value at the boundaries
    #1 reps the second which has restrictions on the derivative at the boundaries
    

    #matrix will be structured so [r1θ1ζ1, r1θ1ζ2, ...r1θ1ζnn, r1θ2ζ1, ..., r1θnmζnn, r2θ1ζ1, .., rnr, θnm, ζnn]

    #but each r point has to have two values, one for each hermite basis.
    #where r_11 is first radial point, first Hermite, r_12 is first radial second Hermite ie the derivative one
    #so it will actually be [rh1θ1ζ1, rh2θ1ζ1, rh1θ1mζ2...]

    segment = (ζind - 1) + Nn * (θind - 1) + Nm * Nn * (rind - 1 + grid_id[hr])

    #the index from 1 is going to change things a bit.
    #can be made clearer but does seem to be working ish.
    return 1 + basis_id[hr] + 2 * segment
    #return 1+2*Nm*Nn*(rind-1) + basis_id[hbind] + 2*Nn*(θind-1) + 2*(ζind-1) + 2*Nm*Nn*grid_id[hbind]

end


"""
    matrix_dim(grids::FSSGridsT)

Computes the dimension of the matrix.
"""
function matrix_dim(grids::FSSGridsT)

    return 2 * grids.r.N * grids.θ.count * grids.ζ.count
end



"""
    index_to_grid(i::Int64, grids::FSSGridsT)

Reverses the transformation from 3d grid to matrix, i.e. maps matrix index to appropriate point in the grid, used for plotting.
"""
function index_to_grid(i::Int64, grids::FSSGridsT)
    r = div(i-1, 2*grids.θ.count*grids.ζ.count) + 1
    θ = mod(div(i-1, 2*grids.ζ.count), grids.θ.count) + 1
    ζ = mod(div(i-1, 2), grids.ζ.count) + 1

    h = mod(i-1, 2) + 1 #ie even or odd

    return r, θ, ζ, h

end


"""
    local_to_global(node::Int64, ξ::Array{Float64}, grid::Array{Float64})

Converts the local grid where the finite elements are defined to the global coordinates.
"""
function local_to_global(node::Int64, ξ::Array{Float64}, grid::Array{Float64})

    #should this be in basis??
    dr = grid[node+1] - grid[node]

    #first we map the (-1, 1) local coords to (0, 1)
    mp = @. (ξ+1) / 2

    #then we sale by the width of grid points
    mp = mp .* dr

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    rglobal = mp .+ grid[node]

    
    return rglobal, dr
end 



"""
    reconstruct_phi(efuncs::Array{ComplexF64, 2}, nevals::Int64, grids::FSSGridsT)

Reconstructs the 1d eigenfunction output back to the 3d grid.
"""
function reconstruct_phi(efuncs::Array{ComplexF64, 2}, nevals::Int64, grids::FSSGridsT)
    phi = zeros(ComplexF64, nevals, grids.r.N, grids.θ.count, grids.ζ.count)
    #maybe one day we will want dphidr???
    #note that this is the same for both!

    for i in 1:2:matrix_dim(grids)

        #note these are the indicies.
        r, θ, ζ, hs = index_to_grid(i, grids)
        phi[:, r, θ, ζ] = efuncs[i, :]

        #if hs == 1
            #may be the wrong way around!
        #    phi[:, r, θ, ζ] = efuncs[i, :]
        #end
    end
    return phi
end



"""
    reconstruct_phi(efuncs::Array{ComplexF64, 2}, nevals::Int64, grids::FSSGridsT)

Reconstructs the 1d eigenfunction output back to the 3d grid. Single efunc case.
"""
function reconstruct_phi(efunc::Array{ComplexF64}, grids::FSSGridsT)
    phi = zeros(ComplexF64, grids.r.N, grids.θ.count, grids.ζ.count)
    #maybe one day we will want dphidr???
    #note that this is the same for both!

    for i in 1:2:matrix_dim(grids)

        #note these are the indicies.
        r, θ, ζ, hs = index_to_grid(i, grids)
        phi[r, θ, ζ] = efunc[i]

        #if hs == 1
            #may be the wrong way around!
        #    phi[:, r, θ, ζ] = efuncs[i, :]
        #end
    end
    return phi
end



"""
Converts the grid point to the appropriate index in the matrix, for the continuum case.

# Args
- These are inconsistent names compared to the normal case!
"""
function cont_grid_to_index(m, n, nn)
    #TODO -> should probably be in its own case.
    #ie matrix goes
    #[m1n1 m1n2...m1nn, m2n1, m2n2...]
    return n + (m-1) * nn

end

"""
Don't think we use this anymroe!
Computes the number of data points computed in our main loop, used for pre-allocating memory.
#not perfectly clear why this is the way it is.
#extremely complicated function for something so simple.
#still a bit unclear!
"""
function compute_length(nr, nm, nn)
    #this is the total time the loops iterate
    total = nm^2 * nn^2 * (nr-1) * 16 #16 from 4x4 Hermite.

    #number of times both left and right are in the boundary.
    #left and right are each in it 2 * pmd.count^2 * tmd.count^2 * 4
    #two for r=0 and r=a, * 4 for the four Hermites for right_bound
    #then we subtract the number of double ups which is
    bounds = 2 * 2 * nm^2 * nn^2 * 4 - 2 * nm^2 * nn^2

    #when the indicies are boundary inds and equal, we apply the boundary conditions
    #so this term must be added back
    equal_bounds = 2 * nm * nn

    return total - bounds + equal_bounds
end


