#If we ever want to make our code extremely flexible, we could make these multiple dispatch on new grid types.
#however, that is unlikely to be useful.


function grid_to_index_linear(xind::Int64, yind::Int64, hx::Int64, hy::Int64, Nx::Int64, Ny::Int64)

    xgrid = xind + hx - 1
    ygrid = yind + hy - 1
    #display(xgrid)

    if xgrid == Nx + 1
        #display("per")
        xgrid = 1
    end

    if ygrid == Ny + 1
        ygrid = 1
    end

    ind = Nx * (xgrid - 1) + (ygrid - 1) 

    return ind + 1

end


function grid_to_index_linear(xind::Int64, h::Int64, N::Int64)

    xgrid = xind + h - 1

    return xgrid 
end

const grid_id = [0, 0, 1, 1]
const basis_id = [0, 1, 0, 1]


#same idea as 3d Hermite, but this is for 1d case.
function grid_to_index_cubic(xind::Int64, h::Int64, N::Int64)
    
    return 1 + 2*(xind-1) + basis_id[h] + 2*grid_id[h]
end

function grid_to_index_cubic(xind::Int64, yind::Int64, hx::Int64, hy::Int64, Nx::Int64, Ny::Int64)
    
    seg = mod(yind-1 + grid_id[hy], Ny) + (xind-1 + grid_id[hx]) * Ny
    return 1 + 4 * seg + basis_id[hy] + 2*basis_id[hx]
end

function grid_to_index_cubic(xind::Int64, yind::Int64, zind::Int64, hx::Int64, hy::Int64, hz::Int64, Nx::Int64, Ny::Int64, Nz::Int64)
    seg = mod(zind-1 + grid_id[hz], Nz) + Nz * mod(yind-1+grid_id[hy], Ny) + Ny * Nz * (xind-1 + grid_id[hx])

    return 1 + 8 * seg + basis_id[hz] + 2 * basis_id[hy] + 4 * basis_id[hx]
end



function index_to_grid_linear(i::Int64, Nx)

end

function index_to_grid_linear(i::Int64, Nx, Ny)

    xind = div(i-1, Ny) + 1
    yind = mod((i-1), Ny) + 1

    return xind, yind

end

function index_to_grid_cubic(i::Int64, Nx, Ny)

    xind = div(i-1, 4*Ny) + 1
    yind = mod(div(i-1, 4), Ny) + 1

    #h = 
    h = mod(i-1, 4) + 1

    return xind, yind, h

end

function index_to_grid_cubic(i::Int64, Nx, Ny, Nz)

    xind = div(i-1, 8*Ny*Nz) + 1
    yind = mod(div(i-1, 8*Nz), Ny) + 1
    zind = mod(div(i-1, 8), Nz) + 1

    h = mod(i-1, 8) + 1

    return xind, yind, zind, h

end
