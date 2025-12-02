"""

Module for the Magnetic and density fields used. 
The magnetic field functions fills out the BFieldT struct based on a q-profile and radial perturbations in the form
    Asin(m0 * θ + n0 * φ).
"""
module Fields


using ..Structures


import Elliptic: ellipke
import FunctionWrappers: FunctionWrapper


export init_fields


include("MagneticField.jl")

export compute_B!


include("qProfiles.jl") 

export quadratic_q, cantori_q, island_q, gae_q, damping_q


include("DensityProfiles.jl")

export uniform_dens, damping_dens, gae_dens


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
