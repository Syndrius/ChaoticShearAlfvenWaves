


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

function inst_grids(grids::ContGridsT)

    return inst_grid(grids.r), inst_grid(grids.θ), inst_grid(grids.ζ)
end

#have this and alternatives to cut the shit.
#I guess this can be like the user interface that we will probably never use.
#could even do some custom printing of the struct.
#not sure if this will ever be used.
function init_grid(type::Symbol, method::Symbol; N::Int64, start::Real=NaN, stop::Real=NaN, sep1::Float64=0.0, sep2::Float64=1.0, frac::Float64=0.0, gp::Int64=4)
    #TODO - if we ever push this publically, this is probably a good idea to fix.
    if type==:radial

        if method == :FEM
            return rfem_grid()
        else
            display("Radial grid must be FEM.")
        end

    elseif type == :angular

        if method == :FEM
            return afem_grid()
        elseif method == :SM
            return asm_grid()
        else
            display("Angular grids must be either FEM or SM.")
        end
    else
        display("Only radial and angular grids are accepted.")
    end

end
