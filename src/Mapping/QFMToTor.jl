
function qfm_spectrum_to_tor(dir_base::String, tor_grids::GridsT, surfs_dir::String)

    mkpath(dir_base*"/tor_map/efuncs")
    mkpath(dir_base*"/tor_map/efuncs_ft")

    prob, qfm_grids, _ = inputs_from_file(dir_base)

    surfs = load_object(surfs_dir)

    evals = evals_from_file(dir_base)

    un_inds = load_object(joinpath(dir_base, "unique_inds.jld2"))

    nevals = length(un_inds)

    sgrid, ϑgrid, φgrid = inst_grids(qfm_grids)

    #think this process is indep of flux vs rad, just depends how the surfaces where generated
    rgrid, θgrid, ζgrid = inst_grids(tor_grids)

    rms = [] #note that we don't know the size of this yet
    tor_ω = []
    mode_labs = Tuple{Int64, Int64}[]


    ϕ_qfm, ϕ_qfmft = PostProcessing.allocate_phi_arrays(qfm_grids, deriv=true)

    ϕ_tor, ϕ_torft = PostProcessing.allocate_phi_arrays(tor_grids, deriv=false)

    plan_qfm = PostProcessing.create_ft_plan(ϕ_qfmft, qfm_grids)
    plan_tor = PostProcessing.create_ft_plan(ϕ_torft, tor_grids)

    rmarray = Array{Int64}(undef, tor_grids.x2.N, tor_grids.x3.N)
    ϕ_tormarray = Array{Float64}(undef, tor_grids.x2.N, tor_grids.x3.N)

    surf_itp, sd = create_surf_itp(surfs)

    CT = CoordTransformT()

    mode_count = 1

    coord_map = qfm_to_tor_coord_map(rgrid, θgrid, ζgrid, CT, surf_itp, sd)

    for i in 1:nevals

        #Not impossible we would only want to map over the chaotic region or something!
        #if evals.x1[i] < sep_min || evals.x1[i] > sep_max
        #    continue
        #end

        #we could generalise this process by creating a get efunc function or something
        #unsure how usefull that would ever be though!
        efunc_read = @sprintf("efunc%05d.hdf5", un_inds[i])
        #unfort doesn't handle complex numbers v well
        efunc_split = load_object(dir_base*"/efuncs_raw/"*efunc_read)

        #ideally this would be preallocated in some way
        efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im

        PostProcessing.reconstruct_phi!(efunc, qfm_grids, ϕ_qfm, ϕ_qfmft, plan_qfm)

        #dont care about maximum in this case
        #amax = argmax(abs.(real.(ϕ_qfm[:, :, 1, 1])))

        
        #bit stupid to pass in the array size here but whatever.
        efunc_map!(ϕ_tor, tor_grids.x1.N, tor_grids.x2.N, tor_grids.x3.N, ϕ_qfm, sgrid, ϑgrid, φgrid, coord_map)
        
        PostProcessing.ft_phi!(ϕ_tor, ϕ_torft, tor_grids, plan_tor)

        rind, mode_lab = PostProcessing.label_mode(ϕ_torft, tor_grids, rmarray, ϕ_tormarray)

        push!(rms, rgrid[rind])
        push!(mode_labs, mode_lab)
        push!(tor_ω, evals.ω[i])

        efunc_write = @sprintf("efunc%05d.jld2", mode_count)

        save_object(dir_base*"/tor_map/efuncs/"*efunc_write, ϕ_tor)
        save_object(dir_base*"/tor_map/efuncs_ft/"*efunc_write, ϕ_torft)
        mode_count += 1
    end

    tor_evals = EvalsT(tor_ω, rms, mode_labs)
    save_object(dir_base*"/tor_map/evals.jld2", tor_evals)
    #so we have a record of the grids used in the mapping
    save_object(dir_base*"/tor_map/grids.jld2", tor_grids)
end


#version without parallel (hopefully!)
#working now!
function qfm_spectrum_to_tor(evals::EvalsT, ϕ_qfm::Array{ComplexF64, 5}, ϕ_qfmft::Array{ComplexF64, 5}, qfm_grids::GridsT, tor_grids::GridsT, surfs::Array{QFMSurfaceT})

    #mkpath(dir_base*"/tor_map/efuncs")
    #mkpath(dir_base*"/tor_map/efuncs_ft")

    #prob, qfm_grids, _ = inputs_from_file(dir_base)

    #surfs = load_object(surfs_dir)

    #evals = evals_from_file(dir_base)

    #un_inds = load_object(joinpath(dir_base, "unique_inds.jld2"))

    nevals = length(evals.ω)

    sgrid, ϑgrid, ζgrid = inst_grids(qfm_grids)

    #think this process is indep of flux vs rad, just depends how the surfaces where generated
    ψgrid, θgrid, φgrid = inst_grids(tor_grids)

    ϕ_tor, ϕ_torft = PostProcessing.allocate_phi_arrays(tor_grids, nevals, deriv=false)

    ϕp, ϕp_ft = PostProcessing.allocate_phi_arrays(tor_grids, deriv=false)

    plan_tor = PostProcessing.create_ft_plan(ϕp_ft, tor_grids)

    surf_itp, sd = create_surf_itp(surfs)

    CT = CoordTransformT()

    ψms = []
    mode_labs = Tuple{Int64, Int64}[]
    tor_ω  = []
    ψmarray = Array{Int64}(undef, tor_grids.x2.N, tor_grids.x3.N)
    ϕ_tormarray = Array{Float64}(undef, tor_grids.x2.N, tor_grids.x3.N)

    coord_map = qfm_to_tor_coord_map(ψgrid, θgrid, φgrid, CT, surf_itp, sd)

    for i in 1:nevals

        
        #bit stupid to pass in the array size here but whatever.
        efunc_map!(ϕp, tor_grids.x1.N, tor_grids.x2.N, tor_grids.x3.N, ϕ_qfm[i, :, :, :, :], sgrid, ϑgrid, ζgrid, coord_map)
        
        PostProcessing.ft_phi!(ϕp, ϕp_ft, tor_grids, plan_tor)

        ψind, mode_lab = PostProcessing.label_mode(ϕp_ft, tor_grids, ψmarray, ϕ_tormarray)

        push!(ψms, ψgrid[ψind])
        push!(mode_labs, mode_lab)
        push!(tor_ω, evals.ω[i])
        ϕ_tor[i, :, :, :] .= ϕp
        ϕ_torft[i, :, :, :] .= ϕp_ft

        ϕp .= 0.0
        ϕp_ft .= 0.0

        #save_object(dir_base*"/tor_map/efuncs/"*efunc_write, ϕ_tor)
        #save_object(dir_base*"/tor_map/efuncs_ft/"*efunc_write, ϕ_torft)
        #mode_count += 1
    end

    tor_evals = EvalsT(tor_ω, ψms, mode_labs)
    #save_object(dir_base*"/tor_map/evals.jld2", tor_evals)
    #so we have a record of the grids used in the mapping
    #save_object(dir_base*"/tor_map/grids.jld2", tor_grids)
    return evals, ϕ_tor, ϕ_torft
end


