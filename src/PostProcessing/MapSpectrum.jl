
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
    #ϕ_tor = zeros(ComplexF64, tor_grids.x1.N, tor_grids.x2.N, tor_grids.x3.N, 8)

    ϕ_tor = zeros(ComplexF64, tor_grids.x1.N, tor_grids.x2.N+1, tor_grids.x3.N+1)

    θgridp = LinRange(0, 2π, tor_grids.x2.N+1)
    ζgridp = LinRange(0, 2π, tor_grids.x3.N+1)

    #not going to store the derivative info of the mapping, maybe we will want to someday.
    #this is probably only a good idea if we assume the island grids are entirely inside or entirely outside.
    #which I think we are doing.
    ϕ_isl, ϕ_islft = allocate_phi_arrays(isl_grids, deriv=false)

    plan = create_ft_plan(ϕ_islft, isl_grids)

    #arrays to store the maximum value of ϕft and the correspondning κ value.
    κmarray = Array{Int64}(undef, isl_grids.x2.N, isl_grids.x3.N)
    ϕ_islmarray = Array{Float64}(undef, isl_grids.x2.N, isl_grids.x3.N)

    #0.0 is the widest part of the island.
    #islands as an array is cooked af here.
    sep_min, sep_max = sepratrix(0.0, prob.isls[1])

    mode_count = 1

    for i in 1:length(evals.ω)

        #this eigenfunction is not an island mode, so we ignore.
        #this is not a good enough check tbh.
        #this just checks that the peak in within the max range of the island 
        #i.e. could peak just outside at the sepratrix.
        #this is just an initial test to remove modes that are fkn ages away!
        if evals.x1[i] < sep_min || evals.x1[i] > sep_max
            continue
        end


        #efunc_read = @sprintf("efunc%05d.hdf5", i)
        #unfort doesn't handle complex numbers v well
        #efunc_split = load_object(dir_base*"/efuncs_raw/"*efunc_read)

        #ideally this would be preallocated in some way
        #efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im
        #raw_dat = @sprintf(dir_base*"/efuncs_raw/efunc%05d.hdf5", i)

        #not sure how best to do this.
        #ϕ_pre .= load_object(raw_dat)

        #this could be odd, as we want to the full thing, but we don't care about the fft.
        #already got it implemented bby
        #reconstruct_phi!(efunc, tor_grids, ϕ_tor)

        efunc_read = @sprintf("efunc%05d.jld2", i)
        ϕ_tor[:, 1:end-1, 1:end-1] .= load_object(dir_base*"/efuncs/"*efunc_read)
        ϕ_tor[:, :, end] = ϕ_tor[:, :, 1]
        ϕ_tor[:, end, :] = ϕ_tor[:, 1, :]
        #is this enough? don't think so!
        ϕ_tor[:, end, end] = ϕ_tor[:, 1, 1]
        #we unfort have to add the periodicity back in, because interpolations.jl is kinda shit.
        amax = argmax(abs.(real.(ϕ_tor[:, :, 1, 1])))

        rmin, rmax = sepratrix(θgrid[amax[2]], prob.isls[1])

        #now a stronger restriction can be placed

        if rmin >= rgrid[amax[1]] || rgrid[amax[1]] >= rmax
            continue
        end
        
        #this is actually a bit of a disaster for islands, for now I think we will just do the inside and ignore the outside
        #so we have to assume the island_grids haev κ <= 1.
        map_tor_to_isl!(ϕ_isl, κgrid, ᾱgrid, τgrid, ϕ_tor, rgrid, θgridp, ζgridp, prob.isls[1])

        #may need to check the plan is not in place or anything stupid.
        #i.e. maybe we write ϕ_isl to file, then fft in place and do the other stuff
        ϕ_islft .= plan * ϕ_isl

        κind, mode_lab = label_mode(ϕ_islft, isl_grids, κmarray, ϕ_islmarray)

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
