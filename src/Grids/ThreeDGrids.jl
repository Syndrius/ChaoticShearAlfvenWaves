struct Grids3dT{G1 <: GridT, G2 <: GridT, G3 <: GridT} <: GridsT
    x1 :: G1
    x2 :: G2
    x3 :: G3
end

function init_grids(x1::GridT, x2::GridT, x3::GridT)

    return Grids3dT(x1, x2, x3)
end

#shorthand for common grids,
const ContGridsT = Grids3dT{ContGridT, SMGridT, SMGridT}
const FSSGridsT = Grids3dT{RadialFEMGridT, SMGridT, SMGridT}
const FFSGridsT = Grids3dT{RadialFEMGridT, PeriodicFEMGridT, SMGridT}
const FFFGridsT = Grids3dT{RadialFEMGridT, PeriodicFEMGridT, PeriodicFEMGridT}

const periodic_types = [:periodic, :angular, :poloidal, :toroidal, :θ, :ζ, :ᾱ, :τ, :ϑ, :φ]
const radial_types = [:radial, :ψ, :r, :s]
const island_types = [:κ, :island]
const spectral_types = [:spectral, :sm, :fourier, :m, :n]


#feel like this doesn't belong here tbh!
function init_grid(type::Symbol, N::Int64; start::Real=0, stop::Real=-1, gp::Int64=4, sep1::Real=0, sep2::Real=0, frac::Float64=0.0, pf::Int64=0, left_bc::Bool=true, incr::Int64=1, f_quad::Int64=3)

    if type in spectral_types
        return init_sm_grid(N, Int(start), incr=incr, f_quad=f_quad)
    end
    
    if type in periodic_types
        if stop == -1
            stop = 2π
        end
        return PeriodicFEMGridT(N, float(start), float(stop), gp, pf)
    end
    if type in island_types
        if stop == -1
            stop = 1.0
        end
        return RadialFEMGridT(N, float(start), float(stop), sep1, sep2, frac, gp, false)
    end
    if type in radial_types
        if stop == -1
            stop = 1.0
        end
        return RadialFEMGridT(N, float(start), float(stop), sep1, sep2, frac, gp, left_bc)
    end
    if type in [:cont, :rc, :ψc, :continuum, :c]
        if stop == -1
            stop = 1.0
        end
        return ContGridT(N, float(start), float(stop))
    end
    display("Grid type not recognised")
end



"""
    inst_grids(grids::GridsT)

Creates the grids used in the computation from the GridT structures.
"""
function inst_grids(grids::Grids3dT)

    return inst_grid(grids.x1), inst_grid(grids.x2), inst_grid(grids.x3)

end

