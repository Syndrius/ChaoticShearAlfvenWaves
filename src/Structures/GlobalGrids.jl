"""
Struct storing the grids in the FSS case.

### Fields
- r::RFEMGridDataT - Struct storing the finite elements radial grid information.
- θ::SMGridDataT - Struct storing the spectral method poloidal grid information.
- ζ::SMGridDataT - Struct storing the spectral method toroidal grid information.
"""
@kwdef struct FSSGridsT <: GridsT
    r :: RFEMGridDataT
    θ :: SMGridDataT
    ζ :: SMGridDataT
end


"""
Struct storing the grids in the FFS case.

### Fields
- r::RFEMGridDataT - Struct storing the finite elements radial grid information.
- θ::AFEMGridDataT - Struct storing the finite element grid.
- ζ::SMGridDataT - Struct storing the spectral method toroidal grid information.
"""
@kwdef struct FFSGridsT <: GridsT
    r :: RFEMGridDataT
    θ :: AFEMGridDataT
    ζ :: SMGridDataT
end


"""
Struct storing the grids in the FFS case.

### Fields
- r::RFEMGridDataT - Struct storing the finite elements radial grid information.
- θ::AEMGridDataT - Struct storing the finite element grid.
- ζ::AFEMGridDataT - Struct storing the finite element grid information.
"""
@kwdef struct FFFGridsT <: GridsT
    r :: RFEMGridDataT
    θ :: AFEMGridDataT
    ζ :: AFEMGridDataT
end


"""
Struct storing the grids in the FFS case.

### Fields
- r::RFEMGridDataT - Struct storing the finite elements radial grid information.
- θ::AEMGridDataT - Struct storing the finite element grid.
- ζ::AFEMGridDataT - Struct storing the finite element grid information.
"""
@kwdef struct ContGridsT <: GridsT
    r :: ContGridDataT
    θ :: SMGridDataT
    ζ :: SMGridDataT
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
        return rfem_grid(N=N, start=start, stop=stop, sep1=sep1, sep2=sep2, frac=frac, gp=gp, left_bc=left_bc)
    end
    #angular finite element grid
    if type==:af
        if isnan(stop)
            stop = 2π
        end
        return afem_grid(N=N, start=start, stop=stop, gp=gp, pf=pf)
    end
    #angular spectral method grid
    if type==:as
        return asm_grid(N=N, start=start, stop=stop, incr=incr, f_quad=f_quad)
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
    init_grids(rgrid::GridDataT, θgrid::GridDataT, ζgrid::GridDataT)

Initialises the grid structure from the three individual grids.
"""
function init_grids(rgrid::GridDataT, θgrid::GridDataT, ζgrid::GridDataT)

    
    if rgrid isa ContGridDataT
        return ContGridsT(rgrid, θgrid, ζgrid)
    end

    if ζgrid isa SMGridDataT

        if θgrid isa SMGridDataT

            return FSSGridsT(rgrid, θgrid, ζgrid)
        else
            return FFSGridsT(rgrid, θgrid, ζgrid)
        end
    elseif ζgrid isa AFEMGridDataT && θgrid isa AFEMGridDataT
        return FFFGridsT(rgrid, θgrid, ζgrid)
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

    return inst_grid(grids.r), inst_grid(grids.θ), inst_grid(grids.ζ)

end


"""
    inst_grids(grids::ContGridsT)

Creates the grids used in the computation from the GridDataT structures.
"""
function inst_grids(grids::ContGridsT)

    return inst_grid(grids.r), inst_grid(grids.θ), inst_grid(grids.ζ)
end

