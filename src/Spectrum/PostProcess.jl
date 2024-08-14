
#file for post processing the data during execution, as we always end up doing it anyway, and this will help the parallel case work better.

#may need more.
#probbaly want to move this into Io.


function post_process(evals::AbstractArray, efuncs::Array{ComplexF64}, grids::FSSGridsT, geo::GeoParamsT)

    ω = geo.R0 .* sqrt.(evals)

    #think the output cont structure is a wee bit cooked.
    #think we should make it a struct.

    #not sure if we want to reconstruct the non-fourier transformed version of this, seems a bit useless.
    ϕft = reconstruct_phi(efuncs, length(evals), grids)


    rgrid, _, mlist, _, _, nlist, _= instantiate_grids(grids)


    rm = zeros(Int64, grids.θ.count, grids.ζ.count)
    ϕm = zeros(Float64, grids.θ.count, grids.ζ.count)

    #mlabs = zeros(Int64, length(evals))
    #nlabs = zeros(Int64, length(evals))
    modelabs = Tuple{Int, Int}[]
    rmode = zeros(Float64, length(evals))
    

    for i in 1:1:length(evals)

        for j in 1:grids.θ.count, k in 1:grids.ζ.count

            rm[j, k] = argmax(abs.(real.(ϕft[i, :, j, k])))
            ϕm[j, k] = abs.(real.(ϕft[i, rm[j, k], j, k]))
        end

        max_mode = argmax(ϕm)

        #not 100% clear what is going on here, but is similar to cka. have to back it in I think.
        


        #modelabs[1, i] = mlist[max_mode[1]]
        #modelabs[2, i] = nlist[max_mode[2]]

        push!(modelabs, (mlist[max_mode[1]], nlist[max_mode[2]]))
        rmode[i] = rgrid[rm[max_mode]]
    end

    evals = EvalsT(ω, rmode, modelabs)

    return evals, nothing, ϕft


end

function post_process(evals::AbstractArray, efuncs::AbstractArray, grids::FFSGridsT, geo::GeoParamsT)


    ω = geo.R0 .* sqrt.(evals)

    #think the output cont structure is a wee bit cooked.
    #think we should make it a struct.

    ϕ = reconstruct_phi(efuncs, length(evals), grids)


    ϕft = zeros(ComplexF64, size(ϕ))

    
    for i in 1:grids.r.N

        for n in 1:grids.ζ.count
        #hopefully the 2d case works as expected
            ϕft[:, i, :, n] = fft(ϕ[:, i, :, n], [2])
        end

    end

    rgrid, _, _, nlist, _= instantiate_grids(grids)

    rm = zeros(Int64, grids.θ.N, grids.ζ.count)
    ϕm = zeros(Float64, grids.θ.N, grids.ζ.count)

    modelabs = Tuple{Int, Int}[]
    rmode = zeros(Float64, length(evals))

    for i in 1:1:length(ω)

        

        for j in 1:grids.θ.N, k in 1:grids.ζ.count

            rm[j, k] = argmax(abs.(real.(ϕft[i, :, j, k])))
            ϕm[j, k] = abs.(real.(ϕft[i, rm[j, k], j, k]))
        end

        max_mode = argmax(ϕm)

        rmode[i] = rgrid[rm[max_mode]]
        #col[i] = (mlist[max_mode[1]], nlist[max_mode[2]])
        #mod is there to correct the mode labels, -1 for julia indexing.
        #this is obvs flawed as we can only resolve Nθ modes, but this should be fine for the neighbouring modes to pf.
        #still not perf, one particular issue is the negative mode numbers.
        #mlab = mod(max_mode[1]-1 + grids.θ.pf, grids.θ.N)
        mlab = mod(max_mode[1]-1, grids.θ.N)
        #think the extra pf here should centre the modes on pf, only really relavant for small Nθ cases.
        if mlab > grids.θ.N/2 + grids.θ.pf
            mlab = mlab - grids.θ.N
        end

        push!(modelabs, (mlab + grids.θ.pf, nlist[max_mode[2]]))
        rmode[i] = rgrid[rm[max_mode]]
    end

    evals = EvalsT(ω, rmode, modelabs)

    return evals, ϕ, ϕft

end

function post_process(evals::Array{ComplexF64}, efuncs::Array{ComplexF64}, grids::FFFGridsT, geo::GeoParamsT)

    #not sure how we want to deal with file writing in this case, I guess just ignore???
    #may want a warning here for -ve values.
    #think for non-parallel case we will just write the entire ϕ to file.
    #This can be done in spectrum to file...
    ω = geo.R0 .* sqrt.(evals)

    #think the output cont structure is a wee bit cooked.
    #think we should make it a struct.

    ϕ = reconstruct_phi(efuncs, length(evals), grids)


    ϕft = zeros(ComplexF64, size(ϕ))

    
    for i in 1:grids.r.N
        #hopefully the 2d case works as expected
        ϕft[:, i, :, :] = fft(ϕ[:, i, :, :], [2, 3])

    end

    rgrid, _, _, = instantiate_grids(grids)

    #arrays for storing the maximum ϕ value and the corresponding radial location.
    rm = zeros(Int64, grids.θ.N, grids.ζ.N)
    ϕm = zeros(Float64, grids.θ.N, grids.ζ.N)
    modelabs = Tuple{Int, Int}[]


    rmode = zeros(Float64, length(evals))
    

    for i in 1:1:length(evals)

        for j in 1:grids.θ.N, k in 1:grids.ζ.N

            rm[j, k] = argmax(abs.(real.(ϕft[i, :, j, k])))
            ϕm[j, k] = abs.(real.(ϕft[i, rm[j, k], j, k]))
        end

        max_mode = argmax(ϕm)

        #not 100% clear what is going on here, but is similar to cka. have to back it in I think.
        mlab = mod(max_mode[1]-1, grids.θ.N)
        nlab = mod(max_mode[2]-1, grids.ζ.N)
        if mlab > grids.θ.N/2
            mlab = mlab - grids.θ.N
        end
        if nlab > grids.ζ.N/2
            nlab = nlab - grids.ζ.N
        end


        push!(modelabs, (mlab + grids.θ.pf, nlab + grids.ζ.pf))
        rmode[i] = rgrid[rm[max_mode]]
    end

    evals = EvalsT(ω, rmode, modelabs)


    return evals, ϕ, ϕft

end