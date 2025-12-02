"""
    tor_spectrum_to_isl(dir_base::String, isl_grids::GridsT)

Maps the spectrum computed in toroidal coordinates, (ψ, θ, φ), to island coordinates, (κ, ᾱ, τ).
This expects the output from CSAWParallel.
"""
function tor_spectrum_to_isl(dir_base::String, isl_grids::GridsT)

    mkpath(dir_base*"/isl_map/efuncs")
    mkpath(dir_base*"/isl_map/efuncs_ft")

    prob, tor_grids, _ = inputs_from_file(dir_base)

    evals = evals_from_file(dir_base)

    un_inds = load_object(joinpath(dir_base, "unique_inds.jld2"))

    nevals = length(un_inds)

    rgrid, θgrid, φgrid = inst_grids(tor_grids)
    κgrid, ᾱgrid, τgrid = inst_grids(isl_grids)

    isl = prob.fields.isls[1]

    κms = [] 
    isl_ω = []
    mode_labs = Tuple{Int64, Int64}[]

    #mostly use fff so ft is not needed, but this keeps the process general.
    ϕ_tor, ϕ_torft = PostProcessing.allocate_phi_arrays(tor_grids, deriv=true)

    ϕ_isl, ϕ_islft = PostProcessing.allocate_phi_arrays(isl_grids, deriv=false)

    plan_tor = PostProcessing.create_ft_plan(ϕ_torft, tor_grids)
    plan_isl = PostProcessing.create_ft_plan(ϕ_islft, isl_grids)

    #arrays to store the maximum value of ϕft and the correspondning κ value.
    κmarray = Array{Int64}(undef, isl_grids.x2.N, isl_grids.x3.N)
    ϕ_islmarray = Array{Float64}(undef, isl_grids.x2.N, isl_grids.x3.N)

    #0.0 is the widest part of the island.
    sep_min, sep_max = separatrix(0.0, isl)

    mode_count = 1

    #maps the coordinates first to efficiently compute the interpolation
    coord_map = tor_to_isl_coord_map(κgrid, ᾱgrid, τgrid, isl)

    for i in 1:nevals

        #this eigenfunction is not an island mode, so we ignore.
        #initial test to remove modes that are certainly not island modes.
        if evals.x1[i] < sep_min || evals.x1[i] > sep_max
            continue
        end

        efunc_read = @sprintf("efunc%05d.hdf5", un_inds[i])
        efunc_split = load_object(dir_base*"/efuncs_raw/"*efunc_read)

        #ideally this would be preallocated in some way
        efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im

        #process the raw efunc to get the full solution for interpolation
        PostProcessing.reconstruct_phi!(efunc, tor_grids, ϕ_tor, ϕ_torft, plan_tor)

        amax = argmax(abs.(real.(ϕ_tor[:, :, :, 1])))

        rmin, rmax = separatrix(isl.m0*θgrid[amax[2]]+isl.n0*φgrid[amax[3]], isl)

        #now a stronger restriction can be placed
        if rmin >= rgrid[amax[1]] || rgrid[amax[1]] >= rmax
            continue
        end

        efunc_map!(ϕ_isl, isl_grids.x1.N, isl_grids.x2.N, isl_grids.x3.N, ϕ_tor, rgrid, θgrid, φgrid, coord_map)

        PostProcessing.ft_phi!(ϕ_isl, ϕ_islft, isl_grids, plan_isl)

        κind, mode_lab = PostProcessing.label_mode(ϕ_islft, isl_grids, κmarray, ϕ_islmarray)

        push!(κms, κgrid[κind])
        push!(mode_labs, mode_lab)
        push!(isl_ω, evals.ω[i])

        efunc_write = @sprintf("efunc%05d.jld2", mode_count)

        save_object(dir_base*"/isl_map/efuncs/"*efunc_write, ϕ_isl)
        save_object(dir_base*"/isl_map/efuncs_ft/"*efunc_write, ϕ_islft)
        mode_count += 1
    end

    isl_evals = EvalsT(isl_ω, κms, mode_labs)
    save_object(dir_base*"/isl_map/evals.jld2", isl_evals)
    #so we have a record of the grids used in the mapping
    save_object(dir_base*"/isl_map/grids.jld2", isl_grids)
end

