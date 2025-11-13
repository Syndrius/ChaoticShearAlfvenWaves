"""
    grid_to_index(x1ind::Int64, x2ind::Int64, x3ind::Int64, hbind::Int64, grids::FSSGridsT)

Converts the index of our 3d grid to the appropriate place in the matrix.
"""
function grid_to_index(x1ind::Int64, x2ind::Int64, x3ind::Int64, hx1::Int64, grids::FSSGridsT)

    Nm = grids.x2.N
    Nn = grids.x3.N

    #id's for the heremite basis
    #tells us which node the basis belongs to, 0 for left, 1 for right of current element
    #these should probably not be redefined all the time...
    
    #tells us which of the basis vectors is being considered here
    #0 reps the first which has restrictions on the value at the boundaries
    #1 reps the second which has restrictions on the derivative at the boundaries
    

    #matrix will be structured so [r1x21x31, r1x21x32, ...r1x21x3nn, r1x22x31, ..., r1x2nmx3nn, r2x21x31, .., rnr, x2nm, x3nn]

    #but each r point has to have two values, one for each hermite basis.
    #where r_11 is first radial point, first Hermite, r_12 is first radial second Hermite ie the derivative one
    #so it will actually be [rh1x21x31, rh2x21x31, rh1x21mx32...]

    segment = (x3ind - 1) + Nn * (x2ind - 1) + Nm * Nn * (x1ind - 1 + grid_id[hx1])

    
    return 1 + basis_id[hx1] + 2 * segment

end



"""
    grid_to_index(x1ind::Int64, x2ind::Int64, x3ind::Int64, hx1::Int64, hx2::Int64, grids::FFSGridsT)

Converts the index of our 3d grid to the appropriate place in the matrix.
"""
function grid_to_index(x1ind::Int64, x2ind::Int64, x3ind::Int64, hx1::Int64, hx2::Int64, grids::FFSGridsT)

    Nx2 = grids.x2.N
    Nn = grids.x3.N
    
    #matrix is structured by segments following the grid.
    #Follows the pattern (0, 0, 0) -> (0, 0, 1) -> ... -> (0, 0, Nn) -> (0, 1, 0) -> (0, Nx2, Nn) -> (1, 0, 0) etc.

    #Each segment contains the 4 basis function combinations at each node. So each segment is made form 4 indicies.
    #rh1x2h1x3, rh1x22x3, rh2x2h1x3, rh2x2h2x3.

    #(subtract 1 due to julia indexing starting from 1)
    # mod is for periodicity, so that if x2ind==Nx2, it is mapped to x2ind = 1.

    #TODO describe this

    #first we find the appropriate segment for our indicies.
    # - Each increment of x3ind moves us one segment.
    # - Each increment of x2ind moves us Nx3 segments.
    # - Each increment of x1ind moves us Nx2 * Nx3 segments.

    #think this is equivalent to earlier work.
    #not sure if this is correct though!
    x2grid = x2ind + grid_id[hx2]
    if x2grid == Nx2+1
        x2grid = 1
    end
    #    display(mod(x2ind-1 + grid_id[hx2], Nx2))
    #end
    #segment = (x3ind-1) + Nn * mod(x2ind-1 + grid_id[hx2], Nx2) + Nx2 * Nn * (x1ind - 1 + grid_id[hx1]) 
    #the -1 doesn't seem to matter inside the mod...
    segment = (x3ind-1) + Nn * (x2grid-1) + Nx2 * Nn * (x1ind - 1 + grid_id[hx1]) 

    #Once we have the segment, we need to find the correct index.
    return 1 + basis_id[hx2] + 2 * basis_id[hx1] + 4 * segment


end



"""
    grid_to_index(x1ind::Int64, x2ind::Int64, x3ind::Int64, hx1::Int64, hx2::Int64, hx3::Int64, grids::FFFGridsT)

Converts the index of our 3d grid to the appropriate place in the matrix.
"""
function grid_to_index(x1ind::Int64, x2ind::Int64, x3ind::Int64, hx1::Int64, hx2::Int64, hx3::Int64, grids::FFFGridsT)

    Nx2 = grids.x2.N
    Nx3 = grids.x3.N
    
    #matrix is structured by segments following the grid.
    #Follows the pattern (0, 0, 0) -> (0, 0, 1) -> ... -> (0, 0, Nn) -> (0, 1, 0) -> (0, Nx2, Nn) -> (1, 0, 0) etc.

    #Each segment contains the 8 basis function combinations at each node. So each segment is made form 8 indicies.
    #rh1x2h1x3h1, rh1x2h1x3h2, rh1x22x3h1, rh1x2h2x3h2, rh2x2h1x3h1, rh2x2h1x3h2, rh2x2h2x3h1, rh2x2h2x3h2.

    #(subtract 1 due to julia indexing starting from 1)
    # mod is for periodicity, so that if x3ind==Nx3, it is mapped to x3ind = 1.


    #first we find the appropriate segment for our indicies.
    # - Each increment of x3ind moves us one segment.
    # - Each increment of x2ind moves us Nx3 segments.
    # - Each increment of x1ind moves us Nx2 * Nx3 segments.
    segment = mod(x3ind-1 + grid_id[hx3], Nx3) + Nx3 * mod(x2ind-1 + grid_id[hx2], Nx2) + Nx2 * Nx3 * (x1ind - 1 + grid_id[hx1]) 

    #Once we have the segment, we need to find the correct index.
    return 1 + 8 * segment + basis_id[hx3] + 2 * basis_id[hx2] + 4 * basis_id[hx1]
end


"""
    grid_to_index(x2ind::Int64, x3ind::Int64, grids::ContGridsT)

Converts the index of the 2d grid into the proper place in the matrix.
"""
function grid_to_index(x2ind::Int64, x3ind::Int64, grids::ContGridsT)

    Nn = grids.x3.N

    return 1 + (x3ind - 1) + (x2ind-1) * Nn

end


"""
    index_to_grid(i::Int64, grids::FSSGridsT)

Reverses the transformation from 3d grid to matrix, i.e. maps matrix index to appropriate point in the grid, used for plotting.
"""
function index_to_grid(i::Int64, grids::FSSGridsT)
    x1 = div(i-1, 2*grids.x2.N*grids.x3.N) + 1
    x2 = mod(div(i-1, 2*grids.x3.N), grids.x2.N) + 1
    x3 = mod(div(i-1, 2), grids.x3.N) + 1

    h = mod(i-1, 2) + 1 #ie even or odd

    return x1, x2, x3, h

end


"""
    index_to_grid(i::Int64, grids::FFSGridsT)

Reverses the transformation from 3d grid to matrix, i.e. maps matrix index to appropriate point in the grid, used for plotting.
"""
function index_to_grid(i::Int64, grids::FFSGridsT)
    x1 = div(i-1, 4*grids.x2.N*grids.x3.N) + 1
    x2 = mod(div(i-1, 4*grids.x3.N), grids.x2.N) + 1
    x3 = mod(div(i-1, 4), grids.x3.N) + 1

    h = mod(i-1, 4) + 1 #ie even or odd

    return x1, x2, x3, h

end




"""
    index_to_grid(i::Int64, grids::FFFGridsT)

Reverses the transformation from 3d grid to matrix, i.e. maps matrix index to appropriate point in the grid, used for plotting.
"""
function index_to_grid(i::Int64, grids::FFFGridsT)
    x1 = div(i-1, 8*grids.x2.N*grids.x3.N) + 1
    x2 = mod(div(i-1, 8*grids.x3.N), grids.x2.N) + 1
    x3 = mod(div(i-1, 8), grids.x3.N) + 1

    h = mod(i-1, 8) + 1 #ie even or odd

    return x1, x2, x3, h

end



#this doesn't really work because evals are returned in order not based on the matrix structure
#this only really works for evals
function index_to_grid(i::Int64, grids::ContGridsT)

    #guess!
    x2 = mod(div(i-1, grids.x3.N), grids.x2.N) + 1
    x3 = mod(i-1, grids.x3.N) + 1

    return x2, x3

end
