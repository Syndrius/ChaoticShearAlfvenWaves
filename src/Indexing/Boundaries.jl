
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

    if grids.r.left_bc == true
        return vcat(left_boundary1, left_boundary2, left_boundary3, left_boundary4, right_boundary1, right_boundary2, right_boundary3, right_boundary4)
    end

    #island case has no left boundaries.
    #return vcat(right_boundary1, right_boundary2, right_boundary3, right_boundary4)
    return vcat(left_boundary1, left_boundary2, left_boundary3, left_boundary4, right_boundary1, right_boundary2, right_boundary3, right_boundary4)

end


"""
    compute_boundary_inds(grids::FFSGridsT)

Computes the indicies of the matrices that correspond to the radial boundaries. Boundaries are ϕ(a) = ϕ(0) = 0. Periodic boundaries are handled elsewhere.
"""
function compute_boundary_inds(grids::FFSGridsT)
    
    Nr = grids.r.N
    Nθ = grids.θ.N
    Nn = grids.ζ.N
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

    if grids.r.left_bc == true
        return  vcat(left_boundary1, left_boundary2, right_boundary1, right_boundary2)
    end

    return vcat(right_boundary1, right_boundary2)
    #return vcat(left_boundary2, right_boundary1, right_boundary2)
    #return vcat(left_boundary1, left_boundary2)

end


"""
    compute_boundary_inds(grids::FSSGridsT)

Computes the indicies of the matrices that correspond to the radial boundaries. Boundaries are ϕ(a) = ϕ(0) = 0. Periodic boundaries are handled elsewhere.
"""
function compute_boundary_inds(grids::FSSGridsT)

    Nr = grids.r.N
    Nm = grids.θ.N
    Nn = grids.ζ.N

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

    if grids.r.left_bc == true
        return vcat(left_boundary, right_boundary)
    end


    #island case doesn't have this boundary.
    return collect(right_boundary)

    
end