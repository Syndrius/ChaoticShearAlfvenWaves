
# What the hek even is this file???
#this may be a placeholder....


@kwdef struct MapGridsT
    Nκ :: Int64
    κmax :: Real
    Nᾱ :: Int64
    Nφ :: Int64
end

function inst_grids(grids::MapGridsT)
    #subject to change.
    return LinRange(0, grids.κmax, grids.Nκ), periodic_grid(grids.Nᾱ, endpoint=true), periodic_grid(grids.Nφ, endpoint=true)
end



#subject to change, ideally this would match the other cases a bit more....
function mode_label(i::Int64, j::Int64, grid::MapGridsT)

    #this doesn't really conform to the other cases. may ned to be changed.

    #this reflects that fft returns modes as
    #[0, 1, ... N/2, -N/2..., -2, -1]
    if i > grid.Nᾱ/2
        mlab = i - grid.Nᾱ - 1
    else
        mlab = i - 1
    end

    if j > grid.Nφ / 2
        nlab = j - grid.Nφ - 1
    else
        nlab = j - 1
    end

    return mlab, nlab
end
