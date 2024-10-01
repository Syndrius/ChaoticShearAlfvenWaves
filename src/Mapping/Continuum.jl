
#can we reconstruct the continuum in island coordinates??? I think not.

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
    isl_ω = evals.ω[0.4 .< evals.r .< 0.6]

    #probably wont actually be this size.
    κms = Array{Float64}(undef, length(isl_ω))
    mode_labs = Tuple{Int, Int}[] 

    κgrid, ᾱgrid, φgrid = inst_grids(islgrids)

    κmarray = Array{Int64}(undef, islgrids.Nᾱ-1, islgrids.Nφ)
    ϕmarray = Array{Float64}(undef, islgrids.Nᾱ-1, islgrids.Nφ)

    #ω_isl = Array{ComplexF64}(undef, length(\isl_))

    #ideally there would be a ft plan in here somewhere lol.
    #there is a disgusting amount of repeated calculation in this function
    #guess we will see if it works at all before worrying about extending lol.
    for i in 1:1:length(isl_ω)

        efunc_write = @sprintf("efunc%04d.jld2", i)

        ind = find_ind(evals, isl_ω[i])


        ϕ = efunc_from_file(dir = dir_base, ind=ind, ft=false);

        #sign shouldn't matter here, as I think we are just doing κ < 1 atm.
        #will be subject to change. Just like the sign part lol.
        ϕ_isl, ϕ_islfft = tor_to_isl(islgrids, ϕ, grids, prob.isl, -1);

        κind, mode_lab = PostProcessing.label_mode(ϕ_islfft, islgrids, κmarray, ϕmarray)


        κms[i] = κgrid[κind]
        push!(mode_labs, mode_lab)

        save_object(dir_base * "/isl_map/efuncs/"*efunc_write, ϕ_isl)
        save_object(dir_base * "/isl_map/efuncs_ft/"*efunc_write, ϕ_islfft)


    end

    isl_evals = EvalsT(isl_ω, κms, mode_labs)
    save_object(dir_base*"/isl_map/evals.jld2", isl_evals)
end