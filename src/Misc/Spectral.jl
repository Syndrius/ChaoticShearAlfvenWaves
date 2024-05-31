
"""
Struct for storing data on the modes considered with the Fourier spectral method.

# Fields
start::Int64 - First mode
count::Int64 - Number of modes
incr::Int64=1 - Gap between modes, defaults to 1.
f_quad::Int64=3 - Multiple for Fourier quadrature, helps reduce aliasing, defulats to 3.
"""
@kwdef struct ModeDataT
    start :: Int64
    count :: Int64
    incr :: Int64 = 1
    f_quad :: Int64 = 3
end


"""
Function to create the relevant grid information for the spectral method. Computes the number of points in the grid, the list of mode numbers and the grid to be used.

# Args
md::ModeDataT - The mode data used to create the grids.
"""
function spectral_grid(md::ModeDataT)

    npoints = md.count * md.f_quad

    mode_list = (md.start:md.incr:md.start + md.incr * md.count)[1:end-1]

    grid = LinRange(0, 2*π / md.incr, npoints+1)[1:end-1]

    return npoints, mode_list, grid

end