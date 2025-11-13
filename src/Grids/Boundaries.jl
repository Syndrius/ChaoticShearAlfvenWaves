#unsure if we want to change this, this is awfully unclear, but I am not sure the other way works.
#if we want a more general case, we will have to change this.

"""
    compute_boundary_inds(grids::FFFGridsT)

Computes the indicies of the matrices that correspond to the radial boundaries. Boundaries are ϕ(a) = ϕ(0) = 0. Periodic boundaries are handled elsewhere.
"""
function compute_boundary_inds(grids::FFFGridsT)

    Nx1 = grids.x1.N
    Nx2 = grids.x2.N
    Nx3 = grids.x3.N
    #note this only gets r boundaries, x2, x3 is assumed periodic and handled elsewhere.
    
    #8 basis functions are produced from the tensor product of 2 basis functions in 3d.
    #Boundary condition is applied to the first basis function in r
    #this occurs in the first 4/8 3d basis functions.
    left_boundary1 = 1:8:8*Nx2 * Nx3
    left_boundary2 = 2:8:8*Nx2 * Nx3
    left_boundary3 = 3:8:8*Nx2 * Nx3
    left_boundary4 = 4:8:8*Nx2 * Nx3
   
    right_boundary1 = 1 + (Nx1 - 1) * 8 * Nx2 * Nx3:8:8*Nx1 * Nx2 * Nx3
    right_boundary2 = 2 + (Nx1 - 1) * 8 * Nx2 * Nx3:8:8*Nx1 * Nx2 * Nx3
    right_boundary3 = 3 + (Nx1 - 1) * 8 * Nx2 * Nx3:8:8*Nx1 * Nx2 * Nx3
    right_boundary4 = 4 + (Nx1 - 1) * 8 * Nx2 * Nx3:8:8*Nx1 * Nx2 * Nx3

    if grids.x1.left_bc == true
        return vcat(left_boundary1, left_boundary2, left_boundary3, left_boundary4, right_boundary1, right_boundary2, right_boundary3, right_boundary4)
    end

    #island case has no left boundaries.
    return vcat(right_boundary1, right_boundary2, right_boundary3, right_boundary4)

end


"""
    compute_boundary_inds(grids::FFSGridsT)

Computes the indicies of the matrices that correspond to the radial boundaries. Boundaries are ϕ(a) = ϕ(0) = 0. Periodic boundaries are handled elsewhere.
"""
function compute_boundary_inds(grids::FFSGridsT)
    
    Nx1 = grids.x1.N
    Nx2 = grids.x2.N
    Nn = grids.x3.N
    #note this only gets r boundaries, x2 is assumed periodic and handles elsewhere.
    
    #4 basis functions are produced from the tensor product of 2 basis functions in 2d.
    #Boundary condition is applied to the first basis function in r
    #this occurs in the first 2/4 2d basis functions.
    left_boundary1 = 1:4:4*Nx2 * Nn
    left_boundary2 = 2:4:4*Nx2 * Nn
   
    right_boundary1 = 1 + (Nx1 - 1) * 4 * Nx2 * Nn:4:4*Nx1 * Nx2 * Nn
    right_boundary2 = 2 + (Nx1 - 1) * 4 * Nx2 * Nn:4:4*Nx1 * Nx2 * Nn

    if grids.x1.left_bc == true
        return  vcat(left_boundary1, left_boundary2, right_boundary1, right_boundary2)
    end

    return vcat(right_boundary1, right_boundary2)

end


"""
    compute_boundary_inds(grids::FSSGridsT)

Computes the indicies of the matrices that correspond to the radial boundaries. Boundaries are ϕ(a) = ϕ(0) = 0. Periodic boundaries are handled elsewhere.
"""
function compute_boundary_inds(grids::FSSGridsT)

    Nx1 = grids.x1.N
    Nm = grids.x2.N
    Nn = grids.x3.N

    #Boundary condition is applied to the first basis function in r
    #this occurs in the first basis functions.
    left_boundary = 1:2:2*Nm*Nn

    #could probably use grid_to_index for this
    right_boundary = 1+(Nx1-1)*2*Nm*Nn:2:Nx1*2*Nm*Nn

    if grids.x1.left_bc == true
        return vcat(left_boundary, right_boundary)
    end

    #island case doesn't have this boundary.
    return collect(right_boundary)
    
end
