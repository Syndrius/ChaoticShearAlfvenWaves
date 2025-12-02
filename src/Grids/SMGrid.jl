"""
Struct for storing data on the modes considered with the Fourier spectral method.

### Fields
- N::Int64 - Number of modes.
- start::Int64 - First mode.
- stop::Int64 - Last mode.
- incr::Int64 - Gap between modes.
- f_quad::Int64 - Multiplier for Fourier quadrature, helps reduce aliasing.
"""
struct SMGridT <: GridT
    N :: Int64
    start :: Int64
    stop :: Int64
    incr :: Int64 
    f_quad :: Int64
end


"""
    function init_sm_grid(N::Int64, start::Int64; incr::Int64=1, f_quad::Int64=3)

Initialises the grid used for the spectral method.
"""
function init_sm_grid(N::Int64, start::Int64; incr::Int64=1, f_quad::Int64=3)

    stop = start + (N-1) * incr

    return SMGridT(N, start, stop, incr, f_quad)
end


"""
    inst_grid(grid::SMGridT)

Instantiates the grid used for computation.
"""
function inst_grid(grid::SMGridT)
    return range(0,  2π / grid.incr, grid.N * grid.f_quad+1)[1:end-1]
end


"""
    mode_list(gd::SMGridT)

Returns the list of modes in the grid.
"""
function mode_list(gd::SMGridT)

    return gd.start:gd.incr:gd.stop
end


"""
    mode_label(i::Int64, grid::SMGridT)

Returns the mode label from an index.
"""
function mode_label(i::Int64, grid::SMGridT)

    mlist = mode_list(grid)

    return mlist[i]
end


"""
    compute_ifft_grid(grid::SMGridT)

Computes the size of the ift grid, used for consistency throughout.
"""
function ifft_size(grid::SMGridT)
    modelist = mode_list(grid)
    #this causes issues with post-processing.
    if length(modelist) > 20 #should probably be set as a constant.
        return length(modelist)
    else 
        return 20
    end
end
