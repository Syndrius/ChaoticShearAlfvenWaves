
"""
Abstract type for each grid.
"""
abstract type GridDataT end


"""
Struct storing the data for finite element grid.

### Fields
- N::Int64 - Number of grid points.
- sep1::Float64 - Left egde of clustered region.
- sep2::Float64 - Right egde of clustered region.
- frac::Float64=0.0 - Fraction of points in clustered region. Set to zero for no clustering.
- gp::Int64=4 - Number of Gaussian quadrature points to use for numerical integration, defaults to 4.
- pf::Int64=0 - Phase factor.
"""
@kwdef struct FEMGridDataT <: GridDataT
    N :: Int64
    sep1 :: Float64 = 0.0
    sep2 :: Float64 = 0.1
    frac :: Float64 = 0.0
    gp :: Int64 = 4
    pf :: Int64 = 0
end


"""
Struct for storing data on the modes considered with the Fourier spectral method.

### Fields
- start::Int64 - First mode.
- count::Int64 - Number of modes.
- incr::Int64=1 - Gap between modes, defaults to 1.
- f_quad::Int64=3 - Multiple for Fourier quadrature, helps reduce aliasing, defaults to 3.
"""
@kwdef struct SMGridDataT <: GridDataT
    start :: Int64
    count :: Int64
    incr :: Int64 = 1
    f_quad :: Int64 = 3
end

"""
Abstract type for the grids.
"""
abstract type GridsT end

"""
Struct storing the grids in the FSS case.

### Fields
- r::FEMGridDataT - Struct storing the finite elements radial grid information.
- θ::SMGridDataT - Struct storing the spectral method poloidal grid information.
- ζ::SMGridDataT - Struct storing the spectral method toroidal grid information.
"""
@kwdef struct FSSGridsT <: GridsT
    r :: FEMGridDataT
    θ :: SMGridDataT
    ζ :: SMGridDataT
end


"""
Struct storing the grids in the FFS case.

### Fields
- r::FEMGridDataT - Struct storing the finite elements radial grid information.
- θ::FEMGridDataT - Struct storing the finite element grid.
- ζ::SMGridDataT - Struct storing the spectral method toroidal grid information.
"""
@kwdef struct FFSGridsT <: GridsT
    r :: FEMGridDataT
    θ :: FEMGridDataT
    ζ :: SMGridDataT
end

"""
Struct storing the grids in the FFS case.

### Fields
- r::FEMGridDataT - Struct storing the finite elements radial grid information.
- θ::FEMGridDataT - Struct storing the finite element grid.
- ζ::FEMGridDataT - Struct storing the finite element grid information.
"""
@kwdef struct FFFGridsT <: GridsT
    r :: FEMGridDataT
    θ :: FEMGridDataT
    ζ :: FEMGridDataT
end


"""
Struct storing the grids for the continuum case.

### Fields
- r::Array{Float64} - Array storing the radial points.
- θ::SMGridDataT - Struct storing the spectral method poloidal grid information.
- ζ::SMGridDataT - Struct storing the spectral method toroidal grid information.
"""
@kwdef struct ContGridsT <: GridsT
    rN :: Int64
    θ :: SMGridDataT
    ζ :: SMGridDataT
end


"""
    init_grids(rgrid::GridDataT, θgrid::GridDataT, ζgrid::GridDataT)

Initialises the grid structure from the three individual grids.
"""
function init_grids(rgrid::GridDataT, θgrid::GridDataT, ζgrid::GridDataT)

    if ζgrid isa SMGridDataT

        if θgrid isa SMGridDataT

            return FSSGridsT(rgrid, θgrid, ζgrid)
        else
            return FFSGridsT(rgrid, θgrid, ζgrid)
        end
    elseif ζgrid isa FEMGridDataT && θgrid isa FEMGridDataT
        return FFFGridsT(rgrid, θgrid, ζgrid)
    end

    display("Not implemented yet...")
    return 0
end


"""
    init_grids(N::Int64, θgrid::SMGridDataT, ζgrid::SMGridDataT)

Initialises the grid structure for the continuum case.
"""
function init_grids(N::Int64, θgrid::SMGridDataT, ζgrid::SMGridDataT)

    return ContGridsT(N, θgrid, ζgrid)
end


"""
    instantiate_grids(grids::FFFGridsT)

Creates the grids used in the computation from the GridDataT structures.
"""
function instantiate_grids(grids::FFFGridsT)

    if grids.r.frac == 0.0
        #rgrid needs to be collected to be consistent with the clustered grid.
        rgrid = collect(range(0.1, 1, grids.r.N))
    else
        rgrid = clustered_grid(grids.r)
    end

    #any island contribution here??
    θgrid = range(0, 2π, grids.θ.N+1)[1:end-1]
    ζgrid = range(0, 2π, grids.ζ.N+1)[1:end-1]

    return rgrid, θgrid, ζgrid

end


"""
    instantiate_grids(grids::FFSGridsT)

Creates the grids used in the computation from the GridDataT structures.
"""
function instantiate_grids(grids::FFSGridsT)

    if grids.r.frac == 0.0
        #rgrid needs to be collected to be consistent with the clustered grid.
        rgrid = collect(range(0, 1, grids.r.N))
    else
        rgrid = clustered_grid(grids.r)
    end

    θgrid = range(0, 2π, grids.θ.N+1)[1:end-1]

    return rgrid, θgrid, sm_grid(grids.ζ)...

end

"""
    instantiate_grids(grids::FSSGridsT)

Creates the grids used in the computation from the GridDataT structures.
"""
function instantiate_grids(grids::FSSGridsT)

    if grids.r.frac == 0.0
        #rgrid needs to be collected to be consistent with the clustered grid.
        rgrid = collect(range(0, 1, grids.r.N))
    else
        rgrid = clustered_grid(grids.r)
    end

    return rgrid, sm_grid(grids.θ)..., sm_grid(grids.ζ)...

end


"""
    instantiate_grids(grids::ContGridsT)

Creates the grids used in the computation from the GridDataT structures.
"""
function instantiate_grids(grids::ContGridsT)

    rgrid = range(0, 1, grids.rN+1)[2:end] #skip zero for continuum case.

    return rgrid, sm_grid(grids.θ)..., sm_grid(grids.ζ)...

end


"""
    sm_grid(gd::SMGridDataT)

Creates the required data for the spectral method.
"""
function sm_grid(gd::SMGridDataT)

    Npoints = gd.count * gd.f_quad

    mode_list = (gd.start:gd.incr:gd.start + gd.incr * gd.count)[1:end-1]

    grid = range(0, 2*π / gd.incr, Npoints+1)[1:end-1]

    return Npoints, mode_list, grid

end


"""
    clustered_grid(grid::FEMGridDataT)

Creates a grid with values clustered between two points for the finite element method.
"""
function clustered_grid(grid::FEMGridDataT)

    #done so old way still works!
    N = grid.N
    sep1 = grid.sep1
    sep2 = grid.sep2
    frac = grid.frac

    nclust = Int(floor(frac*N))

    rclust = LinRange(sep1, sep2, nclust)

    nrest = N-nclust

    nright = Int(floor((1-sep2)*nrest))

    nleft = nrest-nright

    rleft = LinRange(0.1, sep1, nleft+1)[1:end-1]

    rright = LinRange(sep2, 1, nright+1)[2:end]

    return vcat(rleft, rclust, rright)
end


"""
    init_fem_grid(;N::Int64, sep1::Float64=0.0, sep2::Float64=1.0, frac::Float64=0.0, gp::Int64=4, pf::Int64=0)

Creates the grid needed for finite elements.
"""
function init_fem_grid(;N::Int64, sep1::Float64=0.0, sep2::Float64=1.0, frac::Float64=0.0, gp::Int64=4, pf::Int64=0)

    return FEMGridDataT(N, sep1, sep2, frac, gp, pf)
end


"""
    init_sm_grid(;start::Int64, count::Int64, incr::Int64=1, f_quad=3)

Creates the grid needed for the spectral method.
"""
function init_sm_grid(;start::Int64, count::Int64, incr::Int64=1, f_quad=3)

    return SMGridDataT(start, count, incr, f_quad)
end


"""
    mode_label(i::Int64, grid::FEMGridDataT)

Converts a fourier transformed grid point into the proper mode label.
"""
function mode_label(i::Int64, grid::FEMGridDataT)

    #subtract pf here to return the grid to what it would have been 
    #without pf.
    #still a bit sketchy.
    mlab = rem(i-1-grid.pf, grid.N)

    #this reflects that fft returns modes as
    #[0, 1, ... N/2, -N/2..., -2, -1]
    if mlab > grid.N/2
        mlab = mlab - grid.N
    end

    return mlab + grid.pf
end

function mode_label(i::Int64, grid::SMGridDataT)

    _, mlist, _ = sm_grid(grid)

    return mlist[i]
end


#mainly just so we can have consistency for plotting etc.
#should probbaly just return the grid???
#unsure if this is just fake when plotting, but seems to work so oh well.
function compute_ifft_grid(grid::SMGridDataT)
    _, modelist, _ = sm_grid(grid)
    if length(modelist) > 20
        return length(modelist)
    else 
        return 20
    end
end
