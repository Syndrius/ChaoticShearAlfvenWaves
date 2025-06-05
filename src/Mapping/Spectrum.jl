
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

    #perhaps we should check if the efunc_deriv (or whatever) folder exists then we don't reconstruct?

    prob, tor_grids, _ = inputs_from_file(dir=dir_base)

    #display(prob)

    #pretty stupid that we have to do this
    isl = inst_island(prob.isls[1])
    #so the island was not instantiated properly.
    #note that island mode 21 q is hard coded. pretty shite tbh.

    #temp solution, we may need to make our island creation a bit more robust if we are going to keep doing this
    #isl = inst_island(init_island(m0=2, n0=-1, A=0.001, qp=2.0, r0=0.5))

    evals = evals_from_file(dir=dir_base)

    rgrid, θgrid, ζgrid = inst_grids(tor_grids)

    κgrid, ᾱgrid, τgrid = inst_grids(isl_grids)

    κms = [] #note that we don't know the size of this yet
    isl_ω = []
    mode_labs = Tuple{Int64, Int64}[]
    #stores the preprocessed raw data.
    #ideally this can be used, but I think it doesn't work due to the stupid complex split in petsc.
    #ϕ_pre = zeros(8 * tor_grids.x1.N * tor_grids.x2.N * tor_grids.x3.N) 
    #this will be the original eigenfunction used.
    #the hermite version is slow af. May want to come back to this
    #it should be more accurate, but who knows tbh.
    ϕ_tor = zeros(ComplexF64, tor_grids.x1.N, tor_grids.x2.N, tor_grids.x3.N, 8)

    #ϕ_tor = zeros(ComplexF64, tor_grids.x1.N, tor_grids.x2.N, tor_grids.x3.N)

    #θgridp = LinRange(0, 2π, tor_grids.x2.N+1)
    #ζgridp = LinRange(0, 2π, tor_grids.x3.N+1)

    #not going to store the derivative info of the mapping, maybe we will want to someday.
    #this is probably only a good idea if we assume the island grids are entirely inside or entirely outside.
    #which I think we are doing.
    ϕ_isl, ϕ_islft = PostProcessing.allocate_phi_arrays(isl_grids, deriv=false)

    plan = PostProcessing.create_ft_plan(ϕ_islft, isl_grids)

    #arrays to store the maximum value of ϕft and the correspondning κ value.
    κmarray = Array{Int64}(undef, isl_grids.x2.N, isl_grids.x3.N)
    ϕ_islmarray = Array{Float64}(undef, isl_grids.x2.N, isl_grids.x3.N)

    #0.0 is the widest part of the island.
    #islands as an array is cooked af here.
    sep_min, sep_max = sepratrix(0.0, isl)

    mode_count = 1

    #maps the coordinates first to efficiently compute the interpolation
    coord_map = tor_to_isl_coord_map(κgrid, ᾱgrid, τgrid, isl)

    for i in 1:length(evals.ω)
    #for i in 1:3

        #this eigenfunction is not an island mode, so we ignore.
        #this is not a good enough check tbh.
        #this just checks that the peak in within the max range of the island 
        #i.e. could peak just outside at the sepratrix.
        #this is just an initial test to remove modes that are fkn ages away!
        if evals.x1[i] < sep_min || evals.x1[i] > sep_max
            continue
        end


        efunc_read = @sprintf("efunc%05d.hdf5", i)
        #unfort doesn't handle complex numbers v well
        efunc_split = load_object(dir_base*"/efuncs_raw/"*efunc_read)

        #ideally this would be preallocated in some way
        efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im
        #raw_dat = @sprintf(dir_base*"/efuncs_raw/efunc%05d.hdf5", i)

        #not sure how best to do this.
        #ϕ_pre .= load_object(raw_dat)

        #this could be odd, as we want to the full thing, but we don't care about the fft.
        #already got it implemented bby
        PostProcessing.reconstruct_phi!(efunc, tor_grids, ϕ_tor)

        #efunc_read = @sprintf("efunc%05d.jld2", i)
        #ϕ_tor[:, 1:end-1, 1:end-1] .= load_object(dir_base*"/efuncs/"*efunc_read)
        #ϕ_tor[:, :, end] = ϕ_tor[:, :, 1]
        #ϕ_tor[:, end, :] = ϕ_tor[:, 1, :]
        #ϕ_tor .= load_object(dir_base*"efuncs/"*efunc_read)

        #display(ϕ_tor[1:5, 1:5, 1])
        #is this enough? don't think so!
        #ϕ_tor[:, end, end] = ϕ_tor[:, 1, 1]
        #we unfort have to add the periodicity back in, because interpolations.jl is kinda shit.
        #choosing ζ=1 here is probably ok, but not really sure why we did that.
        #amax = argmax(abs.(real.(ϕ_tor[:, :, 1, 1])))
        #ideally we would do this, and actualy compute the sepratrix based on α, 
        #but I think this is fine.
        amax = argmax(abs.(real.(ϕ_tor[:, :, :, 1])))

        rmin, rmax = sepratrix(isl.m0*θgrid[amax[2]]+isl.n0*ζgrid[amax[3]], isl)

        #display((rmin, rmax))
        #display(isl)
        #display((θgridp[amax[2]], rgrid[amax[1]]))

        #now a stronger restriction can be placed

        if rmin >= rgrid[amax[1]] || rgrid[amax[1]] >= rmax
            continue
        end
        
        #this is actually a bit of a disaster for islands, for now I think we will just do the inside and ignore the outside
        #so we have to assume the island_grids haev κ <= 1.
        #map_tor_to_isl!(ϕ_isl, κgrid, ᾱgrid, τgrid, ϕ_tor, rgrid, θgridp, ζgridp, isl)

        efunc_map!(ϕ_isl, isl_grids.x1.N, isl_grids.x2.N, isl_grids.x3.N, ϕ_tor, rgrid, θgrid, ζgrid, coord_map)


        #may need to check the plan is not in place or anything stupid.
        #i.e. maybe we write ϕ_isl to file, then fft in place and do the other stuff
        ϕ_islft .= plan * ϕ_isl

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

    #perhaps we should check if the efunc_deriv (or whatever) folder exists then we don't reconstruct?

    prob, qfm_grids, _ = inputs_from_file(dir=dir_base)

    surfs = load_object(surfs_dir)

    #display(prob)

    #pretty stupid that we have to do this
    #this shouldn't be needed anymore but who knows
    isl = inst_island(prob.isls[1])
    #so the island was not instantiated properly.
    #note that island mode 21 q is hard coded. pretty shite tbh.

    #temp solution, we may need to make our island creation a bit more robust if we are going to keep doing this
    #isl = inst_island(init_island(m0=2, n0=-1, A=0.001, qp=2.0, r0=0.5))

    evals = evals_from_file(dir=dir_base)

    sgrid, ϑgrid, φgrid = inst_grids(qfm_grids)

    κgrid, ᾱgrid, τgrid = inst_grids(isl_grids)

    κms = [] #note that we don't know the size of this yet
    isl_ω = []
    mode_labs = Tuple{Int64, Int64}[]
    #stores the preprocessed raw data.
    #ideally this can be used, but I think it doesn't work due to the stupid complex split in petsc.
    #ϕ_pre = zeros(8 * tor_grids.x1.N * tor_grids.x2.N * tor_grids.x3.N) 
    #this will be the original eigenfunction used.
    #the hermite version is slow af. May want to come back to this
    #it should be more accurate, but who knows tbh.
    ϕ_qfm = zeros(ComplexF64, qfm_grids.x1.N, qfm_grids.x2.N, qfm_grids.x3.N, 8)

    #ϕ_qfm = zeros(ComplexF64, qfm_grids.x1.N, qfm_grids.x2.N+1, qfm_grids.x3.N+1)

    #ϑgridp = LinRange(0, 2π, qfm_grids.x2.N+1)
    #φgridp = LinRange(0, 2π, qfm_grids.x3.N+1)

    #not going to store the derivative info of the mapping, maybe we will want to someday.
    #this is probably only a good idea if we assume the island grids are entirely inside or entirely outside.
    #which I think we are doing.
    ϕ_isl, ϕ_islft = PostProcessing.allocate_phi_arrays(isl_grids, deriv=false)

    plan = PostProcessing.create_ft_plan(ϕ_islft, isl_grids)

    #arrays to store the maximum value of ϕft and the correspondning κ value.
    κmarray = Array{Int64}(undef, isl_grids.x2.N, isl_grids.x3.N)
    ϕ_islmarray = Array{Float64}(undef, isl_grids.x2.N, isl_grids.x3.N)

    #0.0 is the widest part of the island.
    #islands as an array is cooked af here.
    #this probbaly doesn't work paritcularly well with qfm...
    #this first step is probably ok, as this is just to weed out efuncs that are ages away from the island
    #sep_min, sep_max = sepratrix(0.0, isl)

    #this step just allows for the warping of the qfm surfaces, probably not needed.
    #sep_min -= 0.05
    #sep_max += 0.05

    surf_itp, sd = create_surf_itp(surfs)

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

    display(sep_min)
    display(sep_max)

    #display(evals.x1[1:50])


    mode_count = 1

    coord_map = qfm_to_isl_coord_map(κgrid, ᾱgrid, τgrid, isl, CT, surf_itp, sd)
    for i in 1:length(evals.ω)
    #for i in 80:100


        #this eigenfunction is not an island mode, so we ignore.
        #this is not a good enough check tbh.
        #this just checks that the peak in within the max range of the island 
        #i.e. could peak just outside at the sepratrix.
        #this is just an initial test to remove modes that are fkn ages away!
        if evals.x1[i] < sep_min || evals.x1[i] > sep_max
            continue
        end


        efunc_read = @sprintf("efunc%05d.hdf5", i)
        #unfort doesn't handle complex numbers v well
        efunc_split = load_object(dir_base*"/efuncs_raw/"*efunc_read)

        #ideally this would be preallocated in some way
        efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im
        #raw_dat = @sprintf(dir_base*"/efuncs_raw/efunc%05d.hdf5", i)

        #not sure how best to do this.
        #ϕ_pre .= load_object(raw_dat)

        #this could be odd, as we want to the full thing, but we don't care about the fft.
        #already got it implemented bby
        PostProcessing.reconstruct_phi!(efunc, qfm_grids, ϕ_qfm)

        #efunc_read = @sprintf("efunc%05d.jld2", i)
        #ϕ_qfm[:, 1:end-1, 1:end-1] .= load_object(dir_base*"/efuncs/"*efunc_read)
        #ϕ_qfm[:, :, end] = ϕ_qfm[:, :, 1]
        #ϕ_qfm[:, end, :] = ϕ_qfm[:, 1, :]

        #display(ϕ_tor[1:5, 1:5, 1])
        #is this enough? don't think so!
        #ϕ_qfm[:, end, end] = ϕ_qfm[:, 1, 1]
        #we unfort have to add the periodicity back in, because interpolations.jl is kinda shit.
        amax = argmax(abs.(real.(ϕ_qfm[:, :, 1, 1])))

        #rmin, rmax = sepratrix(θgrid[amax[2]], isl)

        #display((rmin, rmax))
        #display(isl)
        #display((θgridp[amax[2]], rgrid[amax[1]]))

        #now a stronger restriction can be placed

        #we can probbaly change the other version to be a bit more like this
        #unsure on this condition, but basically, the smin/smax are functions of ϑ, 
        #so this checks if the s min/max at the peak is less/greater than the sepratrix.
        if smin[amax[2]] >= sgrid[amax[1]] || sgrid[amax[1]] >= smax[amax[2]]
            continue
        end
        
        #this is actually a bit of a disaster for islands, for now I think we will just do the inside and ignore the outside
        #so we have to assume the island_grids haev κ <= 1.
        #Mapping.map_qfm_to_isl!(ϕ_isl, κgrid, ᾱgrid, τgrid, ϕ_qfm, sgrid, ϑgridp, φgridp, isl, CT, surf_itp, sd)

        efunc_map!(ϕ_isl, isl_grids.x1.N, isl_grids.x2.N, isl_grids.x3.N, ϕ_qfm, sgrid, ϑgrid, φgrid, coord_map)
        #may need to check the plan is not in place or anything stupid.
        #i.e. maybe we write ϕ_isl to file, then fft in place and do the other stuff
        ϕ_islft .= plan * ϕ_isl

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

#not sure where this function should be, but not here thats for damn sure.
function sepratrix(α::Float64, isl::IslandT)

    r2diff = sqrt(isl.w^2*(1-sin(isl.m0*α/2)^2))

    return sqrt(-r2diff + isl.r0^2), sqrt(r2diff+isl.r0^2)
end


#maps the sepratrix in toroidal coords to qfm coords
function map_sepratrix(rmin::Array{Float64}, rmax::Array{Float64}, θgrid::AbstractArray{Float64}, isl::IslandT, CT::CoordTransformT, surf_itp::QFM.SurfaceITPT, sd::TempSurfT)


    #should all be the same length
    smin = zeros(length(rmin))
    smax = zeros(length(rmax))
    ϑsep = zeros(length(θgrid))

    for i in 1:length(θgrid)
        #this is probbaly a wildy inefficient way of doing this
        #but it is only done ones so who cares.
        s, ϑ, _ =  Mapping.tor_coords_to_qfm(rmin[i], θgrid[i], 0.0, CT, surf_itp, sd)
        smin[i] = s
        ϑsep[i] = ϑ #will this be the same in both cases? probably...
        s, ϑ, _ = Mapping.tor_coords_to_qfm(rmax[i], θgrid[i], 0.0, CT, surf_itp, sd)
        smax[i] = s
    end


    return smin, smax, ϑsep

end

