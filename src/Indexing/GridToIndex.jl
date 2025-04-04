"""
    grid_to_index(rind::Int64, θind::Int64, ζind::Int64, hbind::Int64, grids::FSSGridsT)

Converts the index of our 3d grid to the appropriate place in the matrix.
"""
function grid_to_index(rind::Int64, θind::Int64, ζind::Int64, hr::Int64, grids::FSSGridsT)

    Nm = grids.θ.N
    Nn = grids.ζ.N

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

    
    return 1 + basis_id[hr] + 2 * segment

end



"""
    grid_to_index(rind::Int64, θind::Int64, ζind::Int64, hr::Int64, hθ::Int64, grids::FFSGridsT)

Converts the index of our 3d grid to the appropriate place in the matrix.
"""
function grid_to_index(rind::Int64, θind::Int64, ζind::Int64, hr::Int64, hθ::Int64, grids::FFSGridsT)

    Nθ = grids.θ.N
    Nn = grids.ζ.N
    
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

    #think this is equivalent to earlier work.
    #not sure if this is correct though!
    θgrid = θind + grid_id[hθ]
    if θgrid == Nθ+1
        θgrid = 1
    end
    #    display(mod(θind-1 + grid_id[hθ], Nθ))
    #end
    #segment = (ζind-1) + Nn * mod(θind-1 + grid_id[hθ], Nθ) + Nθ * Nn * (rind - 1 + grid_id[hr]) 
    #the -1 doesn't seem to matter inside the mod...
    segment = (ζind-1) + Nn * (θgrid-1) + Nθ * Nn * (rind - 1 + grid_id[hr]) 

    #Once we have the segment, we need to find the correct index.
    return 1 + basis_id[hθ] + 2 * basis_id[hr] + 4 * segment


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
    index_to_grid(i::Int64, grids::FSSGridsT)

Reverses the transformation from 3d grid to matrix, i.e. maps matrix index to appropriate point in the grid, used for plotting.
"""
function index_to_grid(i::Int64, grids::FSSGridsT)
    r = div(i-1, 2*grids.θ.N*grids.ζ.N) + 1
    θ = mod(div(i-1, 2*grids.ζ.N), grids.θ.N) + 1
    ζ = mod(div(i-1, 2), grids.ζ.N) + 1

    h = mod(i-1, 2) + 1 #ie even or odd

    return r, θ, ζ, h

end


"""
    index_to_grid(i::Int64, grids::FFSGridsT)

Reverses the transformation from 3d grid to matrix, i.e. maps matrix index to appropriate point in the grid, used for plotting.
"""
function index_to_grid(i::Int64, grids::FFSGridsT)
    r = div(i-1, 4*grids.θ.N*grids.ζ.N) + 1
    θ = mod(div(i-1, 4*grids.ζ.N), grids.θ.N) + 1
    ζ = mod(div(i-1, 4), grids.ζ.N) + 1

    h = mod(i-1, 4) + 1 #ie even or odd

    return r, θ, ζ, h

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
