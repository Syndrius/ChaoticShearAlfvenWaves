"""
Struct storing the data required for the Radial grid using the Finite Element Method.
Specifally, this grid will use Direchlet boundary conditions, and is assumed to have the domain (0, 1).
This grid has the additional option of clustering.

### Fields
- N::Int64 - Number of grid points.
- start :: Real - Starting point of the grid.
- stop :: Real - Stopping point of the grid.
- sep1::Float64 - Left edge of clustered region.
- sep2::Float64 - Right edge of clustered region.
- frac::Float64=0.0 - Fraction of points in clustered region. Set to zero for no clustering.
- gp::Int64=4 - Number of Gaussian quadrature points to use for numerical integration, defaults to 4.
"""
@kwdef struct RadFEMGridDataT <: GridDataT
    N :: Int64
    start :: Real = 0
    stop :: Real = 1
    sep1 :: Float64 = 0.0
    sep2 :: Float64 = 0.0
    frac :: Float64 = 0.0
    gp :: Int64 = 4
    left_bc :: Bool = true #false for island case. Perhaps this is a bad name, perhaps we should also condier all boundary options, even derivs???
end



"""
Struct storing the data required for the Angular grid using the Finite Element Method.
Specifally, this grid will use Periodic boundary conditions, and is assumed to have the domain (0, 2π).
This grid has the additional option of a phase factor (pf), to focus on a specifc mode number.

### Fields
- N::Int64 - Number of grid points.
- start :: Real - Starting point of the grid.
- stop :: Real - Stopping point of the grid.
- gp::Int64=4 - Number of Gaussian quadrature points to use for numerical integration, defaults to 4.
- pf::Int64=0 - Phase factor.
"""
@kwdef struct AngFEMGridDataT <: GridDataT
    N :: Int64
    start :: Real = 0
    stop :: Real = 2π
    gp :: Int64 = 4
    pf :: Int64 = 0
end


"""
    radial_fem_grid(; N::Int64, start::Real=0, stop::Real=1, sep1::Float64=0.0, sep2::Float64=1.0, frac::Float64=0.0, gp::Int64=4, left_bc=true)

Creates the grid used for the radial dimension with finite elements.
"""
function radial_fem_grid(; N::Int64, start::Real=0, stop::Real=1, sep1::Float64=0.0, sep2::Float64=1.0, frac::Float64=0.0, gp::Int64=4, left_bc=true)

    return RadFEMGridDataT(N, start, stop, sep1, sep2, frac, gp, left_bc)
end



"""
    angular_fem_grid(; N::Int64, start::Real=0, stop::Real=2π, gp::Int64=4, pf::Int64=0)

Creates the grid used for angular grids with finite elements.
"""
function angular_fem_grid(; N::Int64, start::Real=0, stop::Real=2π, gp::Int64=4, pf::Int64=0)

    return AngFEMGridDataT(N, start, stop, gp, pf)
end


"""
    inst_grid(grid::AngFEMGridDataT)

Instantiates the grid used for computation.
"""
function inst_grid(grid::AngFEMGridDataT)

    return periodic_grid(grid.N, start=grid.start, stop=grid.stop)
end


"""
    inst_grid(grid::RadFEMGridDataT)

Instantiates the grid used for computation.
"""
function inst_grid(grid::RadFEMGridDataT)

    if grid.frac == 0.0
        return collect(LinRange(grid.start, grid.stop, grid.N))
    else
        return clustered_grid(grid)
    end
end


"""
    clustered_grid(grid::RadFEMGridDataT)

Creates a grid with values clustered between two points for the finite element method.
"""
function clustered_grid(grid::RadFEMGridDataT)

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

    rleft = LinRange(grid.start, sep1, nleft+1)[1:end-1]

    rright = LinRange(sep2, grid.stop, nright+1)[2:end]

    return vcat(rleft, rclust, rright)
end


"""
    periodic_grid(N, stop=2π, endpoint=false)

Creates a periodic grid without the endpoint.
"""
function periodic_grid(N; start=0, stop=2π, endpoint=false)

    if endpoint
        return range(start, stop, N)
    else
        return range(start, stop, N+1)[1:end-1]
    end
end


"""
    periodic_grid(N, stop=2π, endpoint=false)

Creates a periodic grid without the endpoint.
"""
function periodic_grid(grid::AngFEMGridDataT; endpoint=false)

    if endpoint
        return range(grid.start, grid.stop, grid.N)
    else
        return range(grid.start, grid.stop, grid.N+1)[1:end-1]
    end
end


"""
    mode_label(i::Int64, grid::AngFEMGridDataT)

Converts a fourier transformed grid point into the proper mode label.
"""
function mode_label(i::Int64, grid::AngFEMGridDataT)

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



function mode_list(gd::AngFEMGridDataT)
    #TODO
    return NaN
end


