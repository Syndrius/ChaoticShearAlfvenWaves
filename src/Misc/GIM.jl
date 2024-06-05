

#Grids Indexing and Matrices!


"""
Computes the indicies of the matrices that correspond to the boundaries. Boundaries are ϕ(a) = ϕ(0) = 0.
Although ϕ(0) should be ϕ'(0)=0 for m=1? It is also not a true boundary condition, rather a regularisation condition.

# Args
nr::Int64 - Number of radial grid points
nm::Int64 - Number of poloidal modes, i.e. how many poloidal points are in the matrix.
nn::Int64 - Number of toroidal modes, i.e. how many toroidal points are in the matrix.
"""
function compute_boundary_inds(nr::Int64, nm::Int64, nn::Int64, mlist::Array{Int64})
    left_boundary = 1:2:2*nm*nn
    #tried to make the distinction between nθ and nm clearer.
    if 0 in mlist
        zero_ind = argmin(abs.(mlist))
        left_boundary = collect(left_boundary)
        #more complicated regularization condition
        #shift the m=0 conditions to the derivative,
        for i in (length(mlist) - zero_ind)*nn+1:(zero_ind)*nn
            left_boundary[i] += 1
        end
    end

    #could probably use grid_to_index for this
    right_boundary = 1+(nr-1)*2*nm*nn:2:nr*2*nm*nn

    return vcat(left_boundary, right_boundary)

end


"""
Modifies the matrix storing the basis functions, Φ, to reflect the local derivatives being taken.

#Args
Φ::Array{ComplexF64, 3} - 9x4xgp matrix, stores the 9 derivtaives of the potential for each 4 Hermite basis functions at each gauss point.
H::Array{Float64, 2} - Hermite basis functions.
dH::Array{Float64, 2} - Derivtive of Hermite basis functions.
ddH::Array{Float64, 2} - Second derivtive of Hermite basis functions.
m::Int64 - Poloidal mode number, i.e. scale factor for derivative of Fourier basis.
n::Int64 - Toroidal mode number, i.e. scale factor for derivative of Fourier basis.
jac::Float64 - Jacobian of local to global transformation, not to be confused with geometric Jacobian.
"""
function create_local_basis!(Φ::Array{ComplexF64, 3}, H::Array{Float64, 2}, dH::Array{Float64, 2}, ddH::Array{Float64, 2}, m::Int64, n::Int64, jac::Float64)

    #Φ_s
    @. Φ[1, :, :] = dH / jac
    #Φ_θ
    @. Φ[2, :, :] = H * m * 1im
    #Φ_ζ
    @. Φ[3, :, :] = H * n * 1im
    #Φ_ss
    @. Φ[4, :, :] = ddH / jac^2
    #Φ_sθ
    @. Φ[5, :, :] = dH * m * 1im / jac
    #Φ_sζ   
    @. Φ[6, :, :] = dH * n * 1im / jac
    #Φ_θθ
    @. Φ[7, :, :] = H * (-m^2)
    #Φ_θζ
    @. Φ[8, :, :] = H * (-m*n)
    #Φ_ζζ
    @. Φ[9, :, :] = H * (-n^2)

    #test case with no derivative!
    #@. Φ[10, :, :] = H

end


"""
Converts the index of our 3d grid to the appropriate place in the matrix.
Matrix is structure as [r1θ1ζ1, r1θ1ζ2, ...r1θ1ζnn, r1θ2ζ1, ..., r1θnmζnn, r2θ1ζ1, .., rnr, θnm, ζnn]
But each r point has to have two values, one for each hermite basis.
Where r_11 is first radial point, first Hermite, r_12 is first radial second Hermite ie the derivative one
So it will actually be [r11θ1ζ1, r12θ1ζ1, r11θ1mζ2...]

# Args
rind::Int64 - Index of radial grid
θind::Int64 - Index of poloidal grid
ζind::Int64 - Index of toroidal grid
hbind::Int64 - Index of Hermite basis functions
nm::Int64 - Number of poloidal modes, i.e. how many poloidal points are in the matrix.
nn::Int64 - Number of toroidal modes, i.e. how many toroidal points are in the matrix.
"""
function grid_to_index(rind::Int64, θind::Int64, ζind::Int64, hbind::Int64, nm::Int64, nn::Int64)

    #id's for the heremite basis
    #tells us which node the basis belongs to, 0 for left, 1 for right of current element
    #these should probably not be redefined all the time...
    
    #tells us which of the basis vectors is being considered here
    #0 reps the first which has restrictions on the value at the boundaries
    #1 reps the second which has restrictions on the derivative at the boundaries
    

    #matrix will be structured so [r1θ1ζ1, r1θ1ζ2, ...r1θ1ζnn, r1θ2ζ1, ..., r1θnmζnn, r2θ1ζ1, .., rnr, θnm, ζnn]

    #but each r point has to have two values, one for each hermite basis.
    #where r_11 is first radial point, first Hermite, r_12 is first radial second Hermite ie the derivative one
    #so it will actually be [r11θ1ζ1, r12θ1ζ1, r11θ1mζ2...]

    #the index from 1 is going to change things a bit.
    #can be made clearer but does seem to be working ish.
    return 1+2*nm*nn*(rind-1) + basis_id[hbind] + 2*nn*(θind-1) + 2*(ζind-1) + 2*nm*nn*grid_id[hbind]

end


"""
Reverses the transformation from 3d grid to matrix, i.e. maps matrix index to appropriate point in the grid, used for plotting.

# Args
i::Int64 Matrix index
nm::Int64 - Number of poloidal modes, i.e. how many poloidal points are in the matrix.
nn::Int64 - Number of toroidal modes, i.e. how many toroidal points are in the matrix.
"""
function index_to_grid(i::Int64, nm::Int64, nn::Int64)
    s = div(i-1, 2*nm*nn) + 1
    θ = mod(div(i-1, 2*nn), nm) + 1
    ζ = mod(div(i-1, 2), nn) + 1

    h = mod(i-1, 2) + 1 #ie even or odd

    return s, θ, ζ, h
end

"""
Reconstructs the 1d eigenfunction output back to the 3d grid.

# Args
efuncs<:Array{ComplexF64, 2} - Eigenfunctions that are outputed when matrix is solved. #this type doesn't work, not sure how to specify it since the efuncs may be real sometimes? Probably mostly complex.
nevals::Int64 - Number of eigen functions found.
nr::Int64 - Number of radial grid points
nm::Int64 - Number of poloidal modes, i.e. how many poloidal points are in the matrix.
nn::Int64 - Number of toroidal modes, i.e. how many toroidal points are in the matrix.
"""
function reconstruct_phi(efuncs::Array{ComplexF64, 2}, nevals::Int64, nr::Int64, nm::Int64, nn::Int64)
    phi = zeros(ComplexF64, nevals, nr, nm, nn)
    #maybe one day we will want dphidr???

    for i in 1:2*nr*nm*nn

        #note these are the indicies.
        s, θ, ζ, hs = index_to_grid(i, nm, nn)

        if hs == 1
            #may be the wrong way around!
            phi[:, s, θ, ζ] = efuncs[i, :]
        end
    end
    return phi
end


"""
Creates a grid with values clusered between two points.

# Args
N::Int64 - Size of the grid
sep1::Float64 - Location of the start of the clustered region.
sep2::Float64 - Location of the end of the clustered region.
frac::Float64 - Fraction of N that should be in the clustered region.
"""
function clustered_grid(N::Int64, sep1::Float64, sep2::Float64, frac::Float64)

    nclust = Int(floor(frac*N))

    rclust = LinRange(sep1, sep2, nclust)

    nrest = N-nclust

    nright = Int(floor((1-sep2)*nrest))

    nleft = nrest-nright

    rleft = LinRange(0, sep1, nleft+1)[1:end-1]

    rright = LinRange(sep2, 1, nright+1)[2:end]

    return vcat(rleft, rclust, rright)
end



"""
Converts the grid point to the appropriate index in the matrix, for the continuum case.

# Args
- These are inconsistent names compared to the normal case!
"""
function cont_grid_to_index(m, n, nn)
    #ie matrix goes
    #[m1n1 m1n2...m1nn, m2n1, m2n2...]
    return n + (m-1) * nn

end

"""
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