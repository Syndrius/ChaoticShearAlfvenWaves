"""

Module for the fields used. 
This includes the density and the magnetic field, computed from an input q-profile and input magnetic islands.
"""
module Fields


using ..Structures


using Elliptic
using FunctionWrappers
import FunctionWrappers: FunctionWrapper


export init_fields


#good
include("MagneticField.jl")

export compute_B!

#good
include("qProfiles.jl") 

export quadratic_q
export cantori_q
export island_q
export gae_q
export damping_q

#good
include("DensityProfiles.jl")

export uniform_dens
export damping_dens
export gae_dens


#good
include("Island.jl")

export separatrix

"""
    init_fields(type=:ψ; q::Function=quadratic_q, dens::Function=uniform_dens, isls::Array{<:IslandT}=IslandT[])

Function that initiates the fields used.
"""
function init_fields(type=:ψ; q::Function=quadratic_q, dens::Function=uniform_dens, isls::Array{<:IslandT}=IslandT[], isl::IslandT=no_isl)

    if isempty(isls) && isl != no_isl
        isls = [isl]
    end

    if q == island_q
        if length(isls) != 1
            display("Island q profile is only compatible with a single island.")
            return 
        end
    end

    if type in [:ψ, :flux, :f]
        if length(isls) > 1 && !(isls[1] isa FluxIslandT)
            display("Island does not match radial variable")
            return
        end
        if q == island_q
            inst_isl = inst_island(isls[1])
            if isnan(inst_isl.A) || isnan(inst_isl.w)
                display("Please specify m0, n0, ψ0, qp and either A or w.")
                return 
            end
            return FluxFieldsT(q, q, dens, FluxIslandT[inst_isl])
        end
        if isempty(isls)
            isls = FluxIslandT[]
        end

        return FluxFieldsT(q, q, dens, isls)
    end
    if type in [:r, :radial, :rad, :s]
        if length(isls) > 1 && !(isls[1] isa RadialIslandT)
            display("Island does not match radial variable")
            return
        end
        if q == island_q
            inst_isl = inst_island(isls[1])
            if isnan(inst_isl.A) || isnan(inst_isl.w)
                display("Please specify m0, n0, ψ0, qp and either A or w.")
                return 
            end
            return RadialFieldsT(q, q, dens, RadialIslandT[inst_isl])
        end
        if isempty(isls)
            isls = RadialIslandT[]
        end

        return RadialFieldsT(q, q, dens, isls)
    end
    if type in [:κ, :isl, :island]
        if length(isls) > 1 && !(isls[1] isa CoordIslandT)
            display("Island does not match radial variable")
            return
        end
        if q != island_q
            display("q profile has been set to island_q.")
        end
        inst_isl = inst_island(isls[1])
        if isnan(inst_isl.A) || isnan(inst_isl.w)
            display("Please specify m0, n0, ψ0, qp and either A or w.")
            return 
        end
        return IslandFieldsT(island_q, island_q, dens, CoordIslandT[inst_isl])
    end

end

end
