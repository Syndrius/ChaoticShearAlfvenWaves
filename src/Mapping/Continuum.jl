
#can we reconstruct the continuum in island coordinates??? I think not.

#this is pretty fkn slow now, especially when we are not using the kappa structure at all.
function mapped_continuum(dir_base::String, islgrids::MapGridsT)
    #perhaps we should allow using MID's solutions, but this probably requires parallel to ever be useful.

    #going to put stuff into an island thingo separate from other stuff.
    mkpath(dir_base*"/isl_map/efuncs") 
    mkpath(dir_base*"/isl_map/efuncs_ft")

    prob, grids = inputs_from_file(dir=dir_base);

    evals = evals_from_file(dir=dir_base);

    #restrict the eigenvalues to only those of interest.
    #these edges are pretty fkn arbitrary.
    #think ideally we would instead have a criteria to check if the max is actually inside the island or somthing...
    #this should be ok for now.
    #think ideally we would determine which of all the modes are actually island modes. 
    #think this criteria is a good starting point
    #but we should check the original modes, to see if there maximum is inside the sepratrix.
    #may b v slow.
    #test_isl_ω = evals.ω[0.42 .< evals.r .< 0.58]

    #isl_inds = find_island_modes(dir_base)

    rgrid, θgrid, _ = inst_grids(grids)

    #probably wont actually be this size.
    κms = [] #Array{Float64}(undef, length(isl_ω))
    mode_labs = Tuple{Int, Int}[] 

    isl_ω = []

    κgrid, ᾱgrid, φgrid = inst_grids(islgrids)

    κmarray = Array{Int64}(undef, islgrids.Nᾱ, islgrids.Nφ)
    ϕmarray = Array{Float64}(undef, islgrids.Nᾱ, islgrids.Nφ)
    κmarray = Array{Int64}(undef, islgrids.Nᾱ-1, islgrids.Nφ)
    ϕmarray = Array{Float64}(undef, islgrids.Nᾱ-1, islgrids.Nφ)

    #ω_isl = Array{ComplexF64}(undef, length(\isl_))

    mode_count = 1

    #ideally there would be a ft plan in here somewhere lol.
    #there is a disgusting amount of repeated calculation in this function
    #guess we will see if it works at all before worrying about extending lol.

    #this is now no longer working at all lol
    #zero CAP in any case.
    #only hope is that this is because of divergence at r-> 0...
    #cannot be sure though! perhaps need to run another case with the quadratic islands to see if we get 
    #a clear CAP. This is fkn annoying, for now we will have to restrict those with r < 0.3 probably...
    #perhaps this is also causing issues with the (0, 0) mode...
    #testing this with the much smaller case does show that the A when r->0 is at least a big part of the problem...
    for i in 1:1:length(evals.ω)

        #avoid bad behaviour at the axis. (hopefully...)
        if evals.r[i] < 0.3
            continue
        end
               

        ϕ = efunc_from_file(dir = dir_base, ind=i, ft=false);

    
        amax = argmax(abs.(ϕ[:, :, 1]))

        rmin, rmax = sepratrix(θgrid[amax[2]], prob.isl)


        #only include modes that peak inside the sepratrix
        if rmin >= rgrid[amax[1]] || rgrid[amax[1]] >= rmax
            continue
        end


        #sign shouldn't matter here, as I think we are just doing κ < 1 atm.
        #will be subject to change. Just like the sign part lol.
        ϕ_isl, ϕ_islfft = tor_to_isl(islgrids, ϕ, grids, prob.isl, -1);

        κind, mode_lab = PostProcessing.label_mode(ϕ_islfft, islgrids, κmarray, ϕmarray)

        
        push!(κms, κgrid[κind])

        efunc_write = @sprintf("efunc%04d.jld2", mode_count)
        push!(mode_labs, mode_lab)

        push!(isl_ω, evals.ω[i])

        save_object(dir_base * "/isl_map/efuncs/"*efunc_write, ϕ_isl)
        save_object(dir_base * "/isl_map/efuncs_ft/"*efunc_write, ϕ_islfft)
        mode_count += 1



    end

    isl_evals = EvalsT(isl_ω, κms, mode_labs)
    save_object(dir_base*"/isl_map/evals.jld2", isl_evals)
end



#assumes ζ = 0
#returns the rmax and rmin of the sepratrix at a given value.
function sepratrix(θ, isl)

    #θ = α with ζ = 0!
    res = sqrt(isl.w^2*(1 - sin(isl.m0*θ/2)^2))

    return sqrt(-res + isl.r0^2), sqrt(res+isl.r0^2)

end