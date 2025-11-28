"""
    qfm_spectrum_to_isl(dir_base::String, isl_grids::GridsT, surfs_dir::String)

Maps the spectrum computed in QFM coordinates, (s, ϑ, ζ), to island coordinates, (κ, ᾱ, τ).
"""
function qfm_spectrum_to_isl(dir_base::String, isl_grids::GridsT, surfs_dir::String)

    mkpath(dir_base*"/isl_map/efuncs")
    mkpath(dir_base*"/isl_map/efuncs_ft")

    prob, qfm_grids, _ = inputs_from_file(dir_base)

    isl = prob.fields.isls[1]

    un_inds = load_object(joinpath(dir_base, "unique_inds.jld2"))

    nevals = length(un_inds)

    surfs = load_object(surfs_dir)

    evals = evals_from_file(dir_base)

    sgrid, ϑgrid, φgrid = inst_grids(qfm_grids)

    κgrid, ᾱgrid, τgrid = inst_grids(isl_grids)

    κms = [] #note that we don't know the size of this yet
    isl_ω = []
    mode_labs = Tuple{Int64, Int64}[]

    ϕ_qfm, ϕ_qfmft = PostProcessing.allocate_phi_arrays(qfm_grids, deriv=true)

    ϕ_isl, ϕ_islft = PostProcessing.allocate_phi_arrays(isl_grids, deriv=false)

    plan_qfm = PostProcessing.create_ft_plan(ϕ_qfmft, qfm_grids)
    plan_isl = PostProcessing.create_ft_plan(ϕ_islft, isl_grids)

    κmarray = Array{Int64}(undef, isl_grids.x2.N, isl_grids.x3.N)
    ϕ_islmarray = Array{Float64}(undef, isl_grids.x2.N, isl_grids.x3.N)

    surf_itp, sd = create_surf_itp(surfs)

    #this computes the sepratrix in toroidal coordinates
    #so they can be mapped to qfm, for checking for island modes.
    θgrid = LinRange(0, 2π, qfm_grids.x2.N+1)[1:end-1]
    rmin = zeros(qfm_grids.x2.N)
    rmax = zeros(qfm_grids.x2.N)

    for (i, θ) in enumerate(θgrid)
        #assumes ζ=0, probbaly ok
        sep_min, sep_max = separatrix(θ, isl)
        rmin[i] = sep_min
        rmax[i] = sep_max
    end

    CT = CoordTransformT()
    smin, smax, ϑsep = map_separatrix(rmin, rmax, θgrid, isl, CT, surf_itp, sd)

    sep_min = minimum(smin)
    sep_max = maximum(smax)

    mode_count = 1

    coord_map = qfm_to_isl_coord_map(κgrid, ᾱgrid, τgrid, isl, CT, surf_itp, sd)
    for i in 1:nevals

        #first remove evals far from island
        if evals.x1[i] < sep_min || evals.x1[i] > sep_max
            continue
        end

        efunc_read = @sprintf("efunc%05d.hdf5", un_inds[i])
        #doesn't handle complex numbers well
        efunc_split = load_object(dir_base*"/efuncs_raw/"*efunc_read)

        #ideally this would be preallocated in some way
        efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im

        PostProcessing.reconstruct_phi!(efunc, qfm_grids, ϕ_qfm, ϕ_qfmft, plan_qfm)

        #for qfm we just take ζ=0.0 so sepratrix mapping is easier.
        amax = argmax(abs.(real.(ϕ_qfm[:, :, 1, 1])))

        #this checks if the s min/max at the peak is less/greater than the sepratrix.
        if smin[amax[2]] >= sgrid[amax[1]] || sgrid[amax[1]] >= smax[amax[2]]
            continue
        end
        
        efunc_map!(ϕ_isl, isl_grids.x1.N, isl_grids.x2.N, isl_grids.x3.N, ϕ_qfm, sgrid, ϑgrid, φgrid, coord_map)
        
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

