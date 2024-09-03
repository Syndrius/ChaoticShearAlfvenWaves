"""
    function post_process(evals::AbstractArray, efuncs::Array{ComplexF64}, grids::FSSGridsT, geo::GeoParamsT)

Proccesses the eigenvalues and eigenfunctions into useful forms. Returns the an EvalsT struct, containing the normalised eigenvalues and data needed to reconstruct the continuum. Also returns the potential realigned with the 3d grid, with the angular variables not fourier transformed, and with the angular variabels fourier transformed.
"""
function post_process(evals::AbstractArray, efuncs::Array{ComplexF64}, grids::FSSGridsT, geo::GeoParamsT)

    ω = geo.R0 .* sqrt.(evals)

    ϕft = reconstruct_phi(efuncs, length(evals), grids)

    θgrid_size = compute_ifft_grid(grids.θ)
    ζgrid_size = compute_ifft_grid(grids.ζ)

    #ϕ = zeros(ComplexF64, size(ϕft)) 
    ϕ = zeros(ComplexF64, length(evals), grids.r.N, θgrid_size, ζgrid_size)

    rgrid, _, mlist, _, _, nlist, _= instantiate_grids(grids)

    maxm = maximum(abs.(mlist))
    maxn = maximum(abs.(nlist))

    #if maxm < 20
    #    maxm = 20
    #end
    #if maxn < 20
    #    maxn = 20
    #end

    

    #ϕ = zeros(ComplexF64, length(evals), grids.r.N, 2*maxm+1, 2*maxn+1)
    ftmlist = vcat(0:maxm, -maxm:-1)
    
    ftnlist = vcat(0:maxn, -maxn:-1)

    ft_array = zeros(ComplexF64, 2*maxm+1, 2*maxn+1)
    """
    for k in 1:length(evals)
        for l in 1:grids.r.N
            ft_array = zeros(ComplexF64, 2*maxm+1, 2*maxn+1)

            for i in 1:1:length(mlist)
                for j in 1:1:length(nlist)

                    indm = maxm + mlist[i] + 1
                    indn = maxn + nlist[j] + 1

                    ft_array[indm, indn] = ϕft[k, l, i, j]

                end
            end
            ϕ[k, l, :, :] = ifft(ft_array)

        end
    end
    """


    #TODO, cannot just iffft. need to think about what mode labels we actually have
    #think we probably don't actually want to use ifft at all,
    #but we have done this before and loss the phase info,
    #so need to be careful.
    #for i in 1:grids.r.N

        #this works, assuming the mode list is symmetric...
    #    ϕ[:, i, :, :] = ifft(ifftshift(ϕft[:, i, :, :], [2]), [2])

    #end

    θgrid = LinRange(0, 2π, θgrid_size+1)[1:end-1]
    ζgrid = LinRange(0, 2π, ζgrid_size+1)[1:end-1]

    #this works, and for general list of m's
    #probably not as efficeint though!
    #think we are better off creating the overkill array and taking ifft.
    
    for i in 1:length(evals)
        for j in 1:grids.r.N
            #for k in 1:grids.ζ.count
            for k in 1:1:length(nlist)
                for l in 1:1:length(mlist)

                    ϕ[i, j, :, :] += ϕft[i, j, l, k] .* exp.(1im * mlist[l] .* θgrid .+ 1im * nlist[k] .* ζgrid' )
                end
            end
        end
    end
    



    


    rm = zeros(Int64, grids.θ.count, grids.ζ.count)
    ϕm = zeros(Float64, grids.θ.count, grids.ζ.count)

    
    modelabs = Tuple{Int, Int}[]
    rmode = zeros(Float64, length(evals))

    for i in 1:1:length(evals)

        for j in 1:grids.θ.count, k in 1:grids.ζ.count

            rm[j, k] = argmax(abs.(real.(ϕft[i, :, j, k])))
            ϕm[j, k] = abs.(real.(ϕft[i, rm[j, k], j, k]))
        end

        max_mode = argmax(ϕm)

        mlab = mlist[max_mode[1]]
        nlab = nlist[max_mode[2]]

        push!(modelabs, (mlab, nlab))
        rmode[i] = rgrid[rm[max_mode]]
    end

    evals = EvalsT(ω, rmode, modelabs)

    return evals, ϕ, ϕft


end


"""
    function post_process(evals::AbstractArray, efuncs::Array{ComplexF64}, grids::FFSGridsT, geo::GeoParamsT)

Proccesses the eigenvalues and eigenfunctions into useful forms. Returns the an EvalsT struct, containing the normalised eigenvalues and data needed to reconstruct the continuum. Also returns the potential realigned with the 3d grid, with the angular variables not fourier transformed, and with the angular variabels fourier transformed.
"""
function post_process(evals::AbstractArray, efuncs::AbstractArray, grids::FFSGridsT, geo::GeoParamsT)


    ω = geo.R0 .* sqrt.(evals)

    ϕfss = reconstruct_phi(efuncs, length(evals), grids)

    ζgrid_size = compute_ifft_grid(grids.ζ)


    ϕ = zeros(ComplexF64, length(evals), grids.r.N, grids.θ.N, ζgrid_size)
    #ϕ = zeros(ComplexF64, size(ϕ_ffs))

    ϕft = fft(ϕfss, [3])

    ζgrid = LinRange(0, 2π, ζgrid_size+1)[1:end-1]

    rgrid, _, _, nlist, _= instantiate_grids(grids)

    #this is kind of pointless, given we don't have any way to look
    #at the ζ components yet
    #but future proofing I guess.
    #perhaps we should write our own function for this weird custom ifft.
    for i in 1:length(evals)
        for j in 1:grids.r.N
            
            for k in 1:grids.θ.N
                for l in 1:1:length(nlist)

                    ϕ[i, j, k, :] += ϕfss[i, j, k, l] .* exp.(1im * nlist[l] .* ζgrid)
                end
            end
        end
    end
    

    

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
        
        mlab = mode_label(max_mode[1], grids.θ)
        nlab = nlist[max_mode[2]]

        push!(modelabs, (mlab, nlab))
        rmode[i] = rgrid[rm[max_mode]]
    end

    evals = EvalsT(ω, rmode, modelabs)

    return evals, ϕ, ϕft

end

"""
    function post_process(evals::AbstractArray, efuncs::Array{ComplexF64}, grids::FFFGridsT, geo::GeoParamsT)

Proccesses the eigenvalues and eigenfunctions into useful forms. Returns the an EvalsT struct, containing the normalised eigenvalues and data needed to reconstruct the continuum. Also returns the potential realigned with the 3d grid, with the angular variables not fourier transformed, and with the angular variabels fourier transformed.
"""
function post_process(evals::Array{ComplexF64}, efuncs::Array{ComplexF64}, grids::FFFGridsT, geo::GeoParamsT)

    
    ω = geo.R0 .* sqrt.(evals)

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
        mlab = mode_label(max_mode[1], grids.θ)
        nlab = mode_label(max_mode[2], grids.ζ)


        push!(modelabs, (mlab, nlab))
        rmode[i] = rgrid[rm[max_mode]]
    end

    evals = EvalsT(ω, rmode, modelabs)


    return evals, ϕ, ϕft

end