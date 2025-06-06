#these functions are working but are a complete mess

#function for reconstructing the full spectrum in a different coordinate set
#this currenlty is only used for solution obtained in parallel, where the raw data is written.
#as we need the derivative information which is typically excluded in our post-processing
#may add functionality to write extra information as well, but given raw data is written by parallel anyway, seems like a better approach.

#note that GridsT may be a bit misleading, I think this will only work for FFF atm. 
#as in fss case, ϕ is artifically extended, which doesn't match the grids, 
#probably fixable, but fff is our primary use case.
#fss with large grids is probably fine.
#I think this file belongs in Post-Processing.
#bit awkward because it assumes parallal input tbh.
function tor_spectrum_to_isl(dir_base::String, isl_grids::GridsT)

    mkpath(dir_base*"/isl_map/efuncs")
    mkpath(dir_base*"/isl_map/efuncs_ft")

    prob, tor_grids, _ = inputs_from_file(dir=dir_base)

    evals = evals_from_file(dir=dir_base)

    rgrid, θgrid, ζgrid = inst_grids(tor_grids)
    κgrid, ᾱgrid, τgrid = inst_grids(isl_grids)

    κms = [] #note that we don't know the size of this yet
    isl_ω = []
    mode_labs = Tuple{Int64, Int64}[]

    #mostly use fff so ft is not needed, but this keeps the process general.
    ϕ_tor, ϕ_torft = PostProcessing.allocate_phi_arrays(tor_grids, derov=true)

    ϕ_isl, ϕ_islft = PostProcessing.allocate_phi_arrays(isl_grids, deriv=false)

    plan_tor = PostProcessing.create_ft_plan(ϕ_torft, tor_grids)
    plan_isl = PostProcessing.create_ft_plan(ϕ_islft, isl_grids)

    #arrays to store the maximum value of ϕft and the correspondning κ value.
    κmarray = Array{Int64}(undef, isl_grids.x2.N, isl_grids.x3.N)
    ϕ_islmarray = Array{Float64}(undef, isl_grids.x2.N, isl_grids.x3.N)

    #0.0 is the widest part of the island.
    #islands as an array is cooked af here.
    sep_min, sep_max = sepratrix(0.0, prob.isls[1])

    mode_count = 1

    #maps the coordinates first to efficiently compute the interpolation
    coord_map = tor_to_isl_coord_map(κgrid, ᾱgrid, τgrid, isl)

    for i in 1:length(evals.ω)

        #this eigenfunction is not an island mode, so we ignore.
        #initial test to remove modes that are certainly not island modes.
        if evals.x1[i] < sep_min || evals.x1[i] > sep_max
            continue
        end

        efunc_read = @sprintf("efunc%05d.hdf5", i)
        #unfort doesn't handle complex numbers v well
        efunc_split = load_object(dir_base*"/efuncs_raw/"*efunc_read)

        #ideally this would be preallocated in some way
        efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im


        #process the raw efunc to get the full solution for interpolation
        PostProcessing.reconstruct_phi!(efunc, tor_grids, ϕ_tor, ϕ_torft, plan_tor)

        amax = argmax(abs.(real.(ϕ_tor[:, :, :, 1])))

        rmin, rmax = sepratrix(isl.m0*θgrid[amax[2]]+isl.n0*ζgrid[amax[3]], isl)

        #now a stronger restriction can be placed
        if rmin >= rgrid[amax[1]] || rgrid[amax[1]] >= rmax
            continue
        end
        

        efunc_map!(ϕ_isl, isl_grids.x1.N, isl_grids.x2.N, isl_grids.x3.N, ϕ_tor, rgrid, θgrid, ζgrid, coord_map)


        #may need to check the plan is not in place or anything stupid.
        #i.e. maybe we write ϕ_isl to file, then fft in place and do the other stuff
        PostProcessing.ft_phi!(ϕ_is, ϕ_islft, isl_grids, plan_isl)

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


function qfm_spectrum_to_isl(dir_base::String, isl_grids::GridsT, surfs_dir::String)

    mkpath(dir_base*"/isl_map/efuncs")
    mkpath(dir_base*"/isl_map/efuncs_ft")

    prob, qfm_grids, _ = inputs_from_file(dir=dir_base)

    surfs = load_object(surfs_dir)

    evals = evals_from_file(dir=dir_base)

    sgrid, ϑgrid, φgrid = inst_grids(qfm_grids)

    κgrid, ᾱgrid, τgrid = inst_grids(isl_grids)

    κms = [] #note that we don't know the size of this yet
    isl_ω = []
    mode_labs = Tuple{Int64, Int64}[]


    ϕ_qfm, ϕ_qfmft = PostProcessing.allocate_phi_arrays(isl_grids, deriv=true)

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
        sep_min, sep_max = sepratrix(θ, isl)
        rmin[i] = sep_min
        rmax[i] = sep_max
    end

    CT = CoordTransformT()
    smin, smax, ϑsep = map_sepratrix(rmin, rmax, θgrid, isl, CT, surf_itp, sd)

    sep_min = minimum(smin)
    sep_max = maximum(smax)

    mode_count = 1

    coord_map = qfm_to_isl_coord_map(κgrid, ᾱgrid, τgrid, isl, CT, surf_itp, sd)
    for i in 1:length(evals.ω)

        #first remove evals far from island
        if evals.x1[i] < sep_min || evals.x1[i] > sep_max
            continue
        end

        efunc_read = @sprintf("efunc%05d.hdf5", i)
        #unfort doesn't handle complex numbers v well
        efunc_split = load_object(dir_base*"/efuncs_raw/"*efunc_read)

        #ideally this would be preallocated in some way
        efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im

        PostProcessing.reconstruct_phi!(efunc, qfm_grids, ϕ_qfm, ϕ_qfmft, plan_qfm)

        #for qfm we just take ζ=0.0 so sepratrix mapping is easier.
        amax = argmax(abs.(real.(ϕ_qfm[:, :, 1, 1])))

        #we can probbaly change the other version to be a bit more like this
        #unsure on this condition, but basically, the smin/smax are functions of ϑ, 
        #so this checks if the s min/max at the peak is less/greater than the sepratrix.
        #also ignores variance with ζ.
        if smin[amax[2]] >= sgrid[amax[1]] || sgrid[amax[1]] >= smax[amax[2]]
            continue
        end
        
        efunc_map!(ϕ_isl, isl_grids.x1.N, isl_grids.x2.N, isl_grids.x3.N, ϕ_qfm, sgrid, ϑgrid, φgrid, coord_map)
        
        PostProcessing.ft_phi!(ϕ_is, ϕ_islft, isl_grids, plan_isl)

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



