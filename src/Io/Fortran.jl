
#interface with the parallel implementation written in fortran.
#this will allow us to read the fortran results easily.

#will essentially be a combination of par_post_process and some extra fortran stuff for the inputs.

function fortran_process(dir::String)

    #will need to fix this eventually...
    deriv = false
    #will ignore deriv case for now, will want to implement that eventually...
    mkpath(dir*"/efuncs")
    mkpath(dir*"/efuncs_ft")

    prob, grids = fortran_inputs(dir)

    


    #need to actually fill in the island details here. Causing many issues.

    #write the julia version of the prob and grids
    inputs_to_file(prob=prob, grids=grids, dir=dir)

    #for multiples, we are probably just going to reprocess each time.
    #ideally this won't happen very often.

    #multiples do not work for some reason. really good.

    file_list = readdir(dir)
    vals = ComplexF64[]

    counts = 0
    count_list = [0]
        
    for file in file_list
        if file[1:4] == "vals"
            
            open(dir*file, "r") do file
        
                s = readlines(file)
                #don't love the global here...
                #global vals = Array{ComplexF64}(undef, length(s))
                for line in s
                    subs = split(strip(line, ['(', ')', ' ']), ",")
                    push!(vals, parse(Float64, subs[1]) + 1im * parse(Float64, subs[2]) )
                end

                #counts how many runs were made and how many evals in each.
                
                push!(count_list, length(vals) - count_list[end])

                counts += 1

    
            end
            
        end
    end

    #this will need to change with updated solving method
    """open(dir*"vals00.dat", "r") do file
        
        s = readlines(file)
        #don't love the global here...
        global vals = Array{ComplexF64}(undef, length(s))
        for (i, line) in enumerate(s)
            subs = split(strip(line, ['(', ')', ' ']), ",")
            vals[i] = parse(Float64, subs[1]) + 1im * parse(Float64, subs[2]) 
        end

        #display()

        #parse.(ComplexF64, s)
    end"""

    #strs = readdir(dir*"/efuncs_raw00")
    #not sure if a lock hdf5 is always written?
    #sometimes loc sometimes lock...
    #may be best to just delete the appropriate number of evals
    #or delete the loc vars from efuncs_raw if this happens again.
    #these cases may just be cooked though lol.
    #almost certainly don't need this guard anymore...
    #if strs[end][end-2:end] == "loc"

    #    nevals = Int64(length(strs)/2)
    #else
    #    nevals = length(strs)
    #end
    nevals = length(vals)

    ϕp, ϕpft = PostProcessing.allocate_phi_arrays(grids, deriv=deriv)

    rms = Array{Float64}(undef, nevals)

    plan = PostProcessing.create_ft_plan(ϕpft, grids)

    rgrid = inst_grids(grids)[1]

    rmarray = Array{Int64}(undef, grids.θ.N, grids.ζ.N)
    ϕmarray = Array{Float64}(undef, grids.θ.N, grids.ζ.N)

    
    #rmode = zeros(Float64, nconv)
    #mlab = zeros(Int64, nconv)
    #nlab = zeros(Int64, nconv)
    ω = Array{ComplexF64}(undef, nevals)
    mode_labs = Tuple{Int, Int}[] 

    #ok so I think this will work on gadi
    #not on mac as we dont have hdf5.
    
    for (run, sub_nevals) in enumerate(count_list)
        for i in 1:sub_nevals


            #this is fkn slow af.
            #need to stick with jld2 I think.
            #i-1 to match fortran staring at 0 for some reason
            efunc_read = @sprintf("efunc%05d.hdf5", (i-1))
            efunc_write = @sprintf("efunc%05d.jld2", i+count_list[run-1])

            #unfort doesn't handle complex numbers v well
            # - 2, -1 for including count_list = 0 at first index
            #and -1 for indexing stuff.
            load_dir = @sprintf("/efuncs_raw%02d/", run-2)
            # "/efuncs_raw0" * String(run) * "/"
            efunc_split = load_object(dir*load_dir*efunc_read)#[1, :]

            efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im


            PostProcessing.reconstruct_phi!(efunc, grids, ϕp, ϕpft, plan)
            
            rind, mode_lab = PostProcessing.label_mode(ϕpft, grids, rmarray, ϕmarray)


            #if deriv
            #    rm, modelab = MID.Spectrum.label_mode(ϕft[:, :, :, 1], grids)
            #else
            #    rm, modelab = MID.Spectrum.label_mode(ϕft, grids)
            #end

            
            #phi, phi_ft, rm, modelab = process_efunc(efunc, grids)

            push!(mode_labs, mode_lab)

            rms[i] = rgrid[rind]

            #normalise the eigenvalues
            ω[i] = prob.geo.R0 * sqrt(vals[i])
            if deriv
                save_object(dir * "/efuncs_deriv/"*efunc_write, ϕp)
                save_object(dir * "/efuncs_ft_deriv/"*efunc_write, ϕpft)
            else
                save_object(dir * "/efuncs/"*efunc_write, ϕp)
                save_object(dir * "/efuncs_ft/"*efunc_write, ϕpft)
            end
            
        end
    end
    evals = EvalsT(ω, rms, mode_labs)

    save_object(dir*"/evals.jld2", evals)

end


function fortran_inputs(dir::String)

    #creates the normal inputs from the fortran inputs.
    #i.e. makes the prob and grid structs.
    #subject to change, as the current params for fortran is kind of garbage and it would be better if it was closer to julia.

    file = dir * "inputs.txt"

    f = open(file, "r")
        # line_number

    Nr = 1
    rgp = 4
    sep1 = 0.0
    sep2 = 1.0
    frac = 0.0

    Nθ = 1
    θgp = 4
    θpf = 2

    Nζ = 1
    ζgp = 4
    ζpf = 2

    R0 = 10.0
    
    A = 0.0
    m0 = 1
    n0 = 1
    #this needs to change..
    q_prof = 7

    
    # read till end of file
    for line in  readlines(f)
    
        # read a new / next line for every iteration           
        s = split(line, " ")       
        var = s[1]
        
        #line += 1
        #println("$line . $var")

        #rgrid
        if var == "Nr"
            Nr = parse(Int64, s[2])
        elseif var == "rgp"
            rgp = parse(Int64, s[2])
        elseif var == "sep1"
            sep1 = parse(Float64, s[2])
            display(sep1)
        elseif var == "sep2"
            sep2 = parse(Float64, s[2])
        elseif var == "frac"
            frac = parse(Float64, s[2])
            #θgrid
        elseif var == "Nt"
            Nθ = parse(Int64, s[2])
        elseif var == "tgp"
            θgp = parse(Int64, s[2])
        elseif var == "tpf"
            θpf = parse(Int64, s[2])
            #ζgrid
        elseif var == "Nz"
            Nζ = parse(Int64, s[2])
        elseif var == "tgp"
            ζgp = parse(Int64, s[2])
        elseif var == "zpf"
            ζpf = parse(Int64, s[2])
            #geometry
        elseif var == "R0"
            R0 = parse(Float64, s[2])

            #flr #may need to ignore d0 stuff...
            #may not need this stuff for plotting tbh
        elseif var == "beta"
            #beta = parse(Float64, s[2])
            continue
            #island #subject to change, using A or w etc
        elseif var == "A"
            #pretty fkn annoying, this can't handle d..
            A = parse(Float64, replace(s[2], "d"=>"e"))
        elseif var == "m0"
            m0 = parse(Int64, s[2])
        elseif var == "n0"
            n0 = parse(Int64, s[2])

            #equil profiles, #TODO
        elseif var == "q"
            q_prof = s[2]
        end
        #there are a few others like save/read mat, but they probbaly don't matter here.

        #this won't handle cases where some of this doesn't exist, eg non-clustered grid
        
    end

    close(f)
        


    rgrid = rfem_grid(N=Nr, sep1=sep1, sep2=sep2, frac=frac)
    θgrid = afem_grid(N=Nθ, pf=θpf)
    ζgrid = afem_grid(N=Nζ, pf=ζpf)

    grids = init_grids(rgrid, θgrid, ζgrid)

    geo = GeoParamsT(R0=R0)
    #island is annoying if not set up properly
    #isl = init_island(m0=2, n0=-1, A=A)

    #current default is just Axel_q, not ideal...
    #going with default metric atm as well. this is not ideal
    prob = init_problem(q=Axel_q, geo=geo)#, isl=isl)

    return prob, grids



end
