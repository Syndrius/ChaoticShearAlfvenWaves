"""

Module for the equilibrium state of the plasma. This includes the density and the magnetic field, computed from an input q-profile and input magnetic islands.
"""
module Fields


using ..Structures


using Elliptic
using FunctionWrappers
import FunctionWrappers: FunctionWrapper


export init_fields


include("MagneticField.jl")

export BFieldT
export compute_B!
#export compute_B_isl!



include("qProfiles.jl") #once again, this has become a disaster. Ideally, we should just have a few q-profile written, and then perhaps we read a file for extras or something

export q_profile
export fu_dam_q
export cantori_q
export low_shear_q
export low_shear_qfm_q
export qfm_q
export qfm_benchmark_q
export island_q



include("DensityProfiles.jl")

export density_profile
export uniform_dens


#will need to suss if we can just specify q for this?
function init_fields(type=:ψ; q::Function=fu_dam_q, dens::Function=uniform_dens, isls::Array{<:IslandT}=IslandT[])

    if q == island_q
        if length(isls) != 1
            display("Island q profile is only compatible with a single island.")
            return 
        end
    end


    if type in [:ψ, :flux, :f]
        if length(isls) > 1 && !(isls[1] isa FluxIslandT)
            display("Island does not match radial varible")
            return
        end
        if q == island_q
            inst_isl = inst_island(isls[1])
            return FluxFieldsT(q, q, dens, FluxIslandT[inst_isl])
        end
        if isempty(isls)
            isls = FluxIslandT[]
        end

        return FluxFieldsT(q, q, dens, isls)
    end
    if type in [:r, :radial, :rad, :s]
        if length(isls) > 1 && !(isls[1] isa RadialIslandT)
            display("Island does not match radial varible")
            return
        end
        if q == island_q
            inst_isl = inst_island(isls[1])
            return RadialFieldsT(q, q, dens, RadialIslandT[inst_isl])
        end
        if isempty(isls)
            isls = RadialIslandT[]
        end

        return RadialFieldsT(q, q, dens, isls)
    end
    if type in [:κ, :isl, :island]
        if length(isls) > 1 && !(isls[1] isa CoordIslandT)
            display("Island does not match radial varible")
            return
        end
        if q != island_q
            display("q profile has been set to island_q.")
            #inst_isl = inst_island(isls[1])
            #return FluxFieldsT(q, q, dens, FluxIslandT[inst_isl])
        end
        #not sure if this is needed!
        inst_isl = inst_island(isls[1])
        return IslandFieldsT(island_q, island_q, dens, CoordIslandT[inst_isl])
    end

end
end
