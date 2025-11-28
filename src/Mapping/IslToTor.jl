
"""
    isl_spectrum_to_tor(evals::EvalsT, ϕ_isl::Array{ComplexF64, 5}, ϕ_islft::Array{ComplexF64, 5}, isl_grids::GridsT, tor_grids::GridsT, prob::ProblemT)

Converts spectrum computed in island coordinates, (κ, ᾱ, τ), to toroidal coordinates, (ψ, θ, φ).
Expects output from serial computation.
"""
function isl_spectrum_to_tor(evals::EvalsT, ϕ_isl::Array{ComplexF64, 5}, ϕ_islft::Array{ComplexF64, 5}, isl_grids::GridsT, tor_grids::GridsT, prob::ProblemT)

    nevals = length(evals.ω)

    κgrid, ᾱgrid, τgrid = inst_grids(isl_grids)

    ψgrid, θgrid, φgrid = inst_grids(tor_grids)

    ϕ_tor, ϕ_torft = PostProcessing.allocate_phi_arrays(tor_grids, nevals, deriv=false)

    ϕp, ϕp_ft = PostProcessing.allocate_phi_arrays(tor_grids, deriv=false)

    plan_tor = PostProcessing.create_ft_plan(ϕp_ft, tor_grids)

    ψms = []
    mode_labs = Tuple{Int64, Int64}[]
    tor_ω  = []
    ψmarray = Array{Int64}(undef, tor_grids.x2.N, tor_grids.x3.N)
    ϕ_tormarray = Array{Float64}(undef, tor_grids.x2.N, tor_grids.x3.N)

    #pre compute the coordinate map for efficient mapping
    coord_map = isl_to_tor_coord_map(ψgrid, θgrid, φgrid, prob.fields.isls[1])

    for i in 1:nevals

        efunc_map!(ϕp, tor_grids.x1.N, tor_grids.x2.N, tor_grids.x3.N, ϕ_isl[i, :, :, :, :], κgrid, ᾱgrid, τgrid, coord_map)
        
        PostProcessing.ft_phi!(ϕp, ϕp_ft, tor_grids, plan_tor)

        ψind, mode_lab = PostProcessing.label_mode(ϕp_ft, tor_grids, ψmarray, ϕ_tormarray)

        push!(ψms, ψgrid[ψind])
        push!(mode_labs, mode_lab)
        push!(tor_ω, evals.ω[i])
        ϕ_tor[i, :, :, :] .= ϕp
        ϕ_torft[i, :, :, :] .= ϕp_ft

        ϕp .= 0.0
        ϕp_ft .= 0.0

    end

    tor_evals = EvalsT(tor_ω, ψms, mode_labs)
    return tor_evals, ϕ_tor, ϕ_torft
end

