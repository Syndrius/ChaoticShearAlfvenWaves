"""
Struct storing the grids in the FSS case.

### Fields
- x1::RadFEMGridDataT - Struct storing the finite elements radial grid information.
- x2::SMGridDataT - Struct storing the spectral method poloidal grid information.
- x3::SMGridDataT - Struct storing the spectral method toroidal grid information.
"""
@kwdef struct FSSGridsT <: GridsT
    x1 :: RadFEMGridDataT
    x2 :: SMGridDataT
    x3 :: SMGridDataT
end


"""
Struct storing the grids in the FFS case.

### Fields
- x1::RadFEMGridDataT - Struct storing the finite elements radial grid information.
- x2::AngFEMGridDataT - Struct storing the finite element grid.
- x3::SMGridDataT - Struct storing the spectral method toroidal grid information.
"""
@kwdef struct FFSGridsT <: GridsT
    x1 :: RadFEMGridDataT
    x2 :: AngFEMGridDataT
    x3 :: SMGridDataT
end


"""
Struct storing the grids in the FFS case.

### Fields
- x1::RadFEMGridDataT - Struct storing the finite elements radial grid information.
- x2::AEMGridDataT - Struct storing the finite element grid.
- x3::AngFEMGridDataT - Struct storing the finite element grid information.
"""
@kwdef struct FFFGridsT <: GridsT
    x1 :: RadFEMGridDataT
    x2 :: AngFEMGridDataT
    x3 :: AngFEMGridDataT
end


"""
Struct storing the grids in the FFS case.

### Fields
- x1::RadFEMGridDataT - Struct storing the finite elements radial grid information.
- x2::AEMGridDataT - Struct storing the finite element grid.
- x3::AngFEMGridDataT - Struct storing the finite element grid information.
"""
@kwdef struct ContGridsT <: GridsT
    x1 :: ContGridDataT
    x2 :: SMGridDataT
    x3 :: SMGridDataT
end


"""
This may be an improvement, it is always hard to tell.
"""
function init_grid(; type::Symbol=:empty, method::Symbol=:empty, bcs::Symbol=:empty, N::Int64, start::Real=0.0, stop::Real=NaN, sep1::Float64=0.0, sep2::Float64=1.0, frac::Float64=0.0, gp::Int64=4, left_bc=true, pf::Int64=0, incr::Int64=1, f_quad::Int64=3)
    #shorthands for the most typical cases
    #radial finite element grid
    if type==:rf
        if isnan(stop)
            stop = 1.0
        end
        return radial_fem_grid(N=N, start=start, stop=stop, sep1=sep1, sep2=sep2, frac=frac, gp=gp, left_bc=left_bc)
    end
    #angular finite element grid
    if type==:af
        if isnan(stop)
            stop = 2π
        end
        return angular_fem_grid(N=N, start=start, stop=stop, gp=gp, pf=pf)
    end
    #angular spectral method grid
    if type==:as
        return angular_sm_grid(N=N, start=start, stop=stop, incr=incr, f_quad=f_quad)
    end

    if type==:rc
        if isnan(stop)
            stop = 1.0
        end
        return ContGridDataT(N=N, start=start, stop=stop)
    end

    if method == :FEM
        if bcs == :dirichlet
            return rfem_grid()
        elseif bcs == :periodic
            return afem_grid()
        else
            display("Only dirichlet and periodic boundary conditions are supported.")
        end
    elseif method == :SM

        if bcs == :periodic
            return asm_grid()
        else
            display("Only periodic boundary conditions are supported.")
        end
    else
        display("Method must be either finite elements (:FEM) or fourier spectral method (:SM).")
    end

end



"""
    init_grids(x1grid::GridDataT, x2grid::GridDataT, x3grid::GridDataT)

Initialises the grid structure from the three individual grids.
"""
function init_grids(x1grid::GridDataT, x2grid::GridDataT, x3grid::GridDataT)

    
    if x1grid isa ContGridDataT
        return ContGridsT(x1grid, x2grid, x3grid)
    end

    if x3grid isa SMGridDataT

        if x2grid isa SMGridDataT

            return FSSGridsT(x1grid, x2grid, x3grid)
        else
            return FFSGridsT(x1grid, x2grid, x3grid)
        end
    elseif x3grid isa AngFEMGridDataT && x2grid isa AngFEMGridDataT
        return FFFGridsT(x1grid, x2grid, x3grid)
    end

    #perhaps make a new radial grid for continuum...
    display("Not implemented yet...")
    return 0
end



"""
    inst_grids(grids::GridsT)

Creates the grids used in the computation from the GridDataT structures.
"""
function inst_grids(grids::GridsT)

    return inst_grid(grids.x1), inst_grid(grids.x2), inst_grid(grids.x3)

end


"""
    inst_grids(grids::ContGridsT)

Creates the grids used in the computation from the GridDataT structures.
"""
function inst_grids(grids::ContGridsT)

    return inst_grid(grids.x1), inst_grid(grids.x2), inst_grid(grids.x3)
end

