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
- frac::Float64 - Fraction of points in clustered region. Set to zero for no clustering.
- gp::Int64 - Number of Gaussian quadrature points to use for numerical integration.
"""
struct RadialFEMGridT <: GridT
    N :: Int64
    start :: Float64
    stop :: Float64
    sep1 :: Float64
    sep2 :: Float64
    frac :: Float64
    gp :: Int64
    left_bc :: Bool
end



"""
Struct storing the data required for the periodic grid using the Finite Element Method.
Specifally, this grid will use Periodic boundary conditions, and is assumed to have the domain (0, 2π).
This grid has the additional option of a phase factor (pf), to focus on a specifc mode number.

### Fields
- N::Int64 - Number of grid points.
- start :: Real - Starting point of the grid.
- stop :: Real - Stopping point of the grid.
- gp::Int64 - Number of Gaussian quadrature points to use for numerical integration, defaults to 4.
- pf::Int64 - Phase factor.
"""
struct PeriodicFEMGridT <: GridT
    N :: Int64
    start :: Float64
    stop :: Float64
    gp :: Int64
    pf :: Int64
end


"""
    init_fem_grid(type::Symbol, N::Int64; start::Real=0.0, stop::Real=-1, bp::Int64=4, sep1::Real=0, sep2::Real=0, frac::Float64=0.0, pf::Int64=0, left_bc::Bool=true)

Initialises a 1d finite element grid.
Can choose a radial (:ψ), periodic (:θ or :φ) or island (:κ) to fill in defualt values.
"""
function init_fem_grid(type::Symbol, N::Int64; start::Real=0.0, stop::Real=-1, bp::Int64=4, sep1::Real=0, sep2::Real=0, frac::Float64=0.0, pf::Int64=0, left_bc::Bool=true)

    if type in periodic_types
        if stop == -1
            stop = 2π
        end
        return PeriodicFEMGridT(N, float(start), float(stop), gp, pf)
    end
    if type in island_types
        if stop == -1
            stop = 0.999 #prevents problems at the sepratrix
        end
        return RadialFEMGridT(N, float(start), float(stop), sep1, sep2, frac, gp, false)
    end
    if type in radial_types
        if stop == -1
            stop = 1.0
        end
        return RadialFEMGridT(N, float(start), float(stop), sep1, sep2, frac, gp, left_bc)
    end
    display("Grid type not recognised")
end



"""
    inst_grid(grid::PeriodicFEMGridT)

Instantiates the grid used for computation.
"""
function inst_grid(grid::PeriodicFEMGridT)

    return periodic_grid(grid.N, start=grid.start, stop=grid.stop)
end


"""
    inst_grid(grid::RadialFEMGridT)

Instantiates the grid used for computation.
"""
function inst_grid(grid::RadialFEMGridT)

    if grid.frac == 0.0
        return collect(LinRange(grid.start, grid.stop, grid.N))
    else
        return clustered_grid(grid)
    end
end


"""
    clustered_grid(grid::RadialFEMGridT)

Creates a grid with values clustered between two points for the finite element method.
"""
function clustered_grid(grid::RadialFEMGridT)

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
function periodic_grid(N::Int64; start=0, stop=2π, endpoint=false)

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
function periodic_grid(grid::PeriodicFEMGridT; endpoint=false)

    if endpoint
        return range(grid.start, grid.stop, grid.N)
    else
        return range(grid.start, grid.stop, grid.N+1)[1:end-1]
    end
end


"""
    mode_label(i::Int64, grid::PeriodicFEMGridT)

Converts a fourier transformed grid point into the proper mode label.
"""
function mode_label(i::Int64, grid::PeriodicFEMGridT)

    #subtract pf here to return the grid to what it would have been 
    #without pf.
    mlab = rem(i-1-grid.pf, grid.N)

    #this reflects that fft returns modes as
    #[0, 1, ... N/2, -N/2..., -2, -1]
    if mlab > grid.N/2
        mlab = mlab - grid.N
    end

    return mlab + grid.pf
end


"""
    mode_list(grid::PeriodicFEMGridT)

Returns the possible modes for a FEM grid.
"""
function mode_list(grid::PeriodicFEMGridT)

    if isodd(grid.N)
        return vcat(0:grid.N ÷ 2,  -grid.N÷2:-1) .+ grid.pf
    else
        return vcat(0:grid.N ÷ 2,  -grid.N÷2+1:-1) .+ grid.pf
    end
end


