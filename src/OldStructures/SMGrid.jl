"""
Struct for storing data on the modes considered with the Fourier spectral method.

### Fields
- N::Int64 - Number of modes.
- start::Int64 - First mode.
- stop::Int64 - Last mode.
- incr::Int64=1 - Gap between modes, defaults to 1.
- f_quad::Int64=3 - Multiple for Fourier quadrature, helps reduce aliasing, defaults to 3.
"""
@kwdef struct SMGridDataT <: GridDataT
    N :: Int64
    start :: Int64
    stop :: Int64
    incr :: Int64 = 1
    f_quad :: Int64 = 3
end


"""
    angular_sm_grid(; N::Int64, start::Int64, stop::Real=NaN, incr::Int64=1, f_quad::Int64=3)

Creates the grid needed for the spectral method.
"""
function angular_sm_grid(; N::Int64, start::Int64, stop::Real=NaN, incr::Int64=1, f_quad::Int64=3)

    #perhaps we should be able to define a start and a stop and N is computed?
    if isnan(stop)
        stop = start + (N - 1) * incr
    end

    return SMGridDataT(N , start, stop, incr, f_quad)
end


"""
    inst_grid(grid::SMGridDataT)

Instantiates the grid used for computation.
"""
function inst_grid(grid::SMGridDataT)

    return range(0,  2π / grid.incr, grid.N * grid.f_quad+1)[1:end-1]

end


"""
    mode_list(gd::SMGridDataT)

Returns the list of modes in the grid.
"""
function mode_list(gd::SMGridDataT)

    return gd.start:gd.incr:gd.stop
end



"""
    mode_label(i::Int64, grid::SMGridDataT)

Returns the mode label from an index.
"""
function mode_label(i::Int64, grid::SMGridDataT)

    mlist = mode_list(grid)

    return mlist[i]
end


"""
    compute_ifft_grid(grid::SMGridDataT)

Computes the size of the ift grid, used for consistency throughout.
#this causes issues for post-processing. Think we can bin this.
"""
function ifft_size(grid::SMGridDataT)
    modelist = mode_list(grid)
    if length(modelist) > 20
        return length(modelist)
    else 
        return 20
    end
end
