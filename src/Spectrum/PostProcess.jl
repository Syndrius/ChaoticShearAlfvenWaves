
#see if we can make post processing better
#and work for both mid and midparallel

#fff is working now.
#still need to consider the derivative case.
function post_process(evals::AbstractArray, efuncs::Array{ComplexF64}, grids::GridsT, geo::GeoParamsT, deriv=false::Bool)

    #deriv is only really valid for fff for now, and probably for a while lol.
    if deriv
        ϕ_g = zeros(ComplexF64, length(evals), phi_size(grids, true)..., 8)
        ϕft_g = zeros(ComplexF64, length(evals), phi_size(grids, false)..., 8)
    else
        ϕ_g = zeros(ComplexF64, length(evals), phi_size(grids, true)...)
        ϕft_g = zeros(ComplexF64, length(evals), phi_size(grids, false)...)
    end

    rms = zeros(length(evals))

    mode_labs = Tuple{Int, Int}[] 

    

    #ϕ_recon = reconstruct_phi(efuncs, length(evals), grids)
    for i in 1:length(evals)

        if deriv
            ϕ_recon = reconstruct_phi_deriv(efuncs[:, i], grids)
        else
            ϕ_recon = reconstruct_phi(efuncs[:, i], grids)
        end

        

        #ϕ, ϕft = ft_phi(ϕ_recon[i, :, :, :], grids)
        ϕ, ϕft = ft_phi(ϕ_recon, grids)

        

        #this is crap.
        #really shouldn't need an extra if condition here.
        if deriv 
            rm, mode_lab = label_mode(ϕft[:, :, :, 1], grids)
            ϕ_g[i, :, :, :, :] = ϕ
            ϕft_g[i, :, :, :, :] = ϕft
        else
            rm, mode_lab = label_mode(ϕft, grids)
            ϕ_g[i, :, :, :] = ϕ
            ϕft_g[i, :, :, :] = ϕft
        end

        rms[i] = rm
        push!(mode_labs, mode_lab)


    end

    ω = geo.R0 .* sqrt.(evals)

    evals = EvalsT(ω, rms, mode_labs)

    return evals, ϕ_g, ϕft_g

end


function ft_phi(ϕ_recon, grids::FSSGridsT)

    θgrid_size = compute_ifft_grid(grids.θ)
    ζgrid_size = compute_ifft_grid(grids.ζ)

    θgrid = LinRange(0, 2π, θgrid_size+1)[1:end-1]
    ζgrid = LinRange(0, 2π, ζgrid_size+1)[1:end-1]

    ϕ = zeros(ComplexF64, grids.r.N, θgrid_size, ζgrid_size)

    _, _, mlist, _, _, nlist, _= instantiate_grids(grids)
    
    for i in 1:grids.r.N
        #for k in 1:grids.ζ.count
        for j in 1:1:length(mlist)
            for k in 1:1:length(nlist)

                ϕ[i, :, :] += ϕ_recon[i, j, k] .* exp.(1im * mlist[j] .* θgrid .+ 1im * nlist[k] .* ζgrid' )
            end
        end
    end

    return ϕ, ϕ_recon
end


function ft_phi(ϕ_recon, grids::FFSGridsT)

    ζgrid_size = compute_ifft_grid(grids.ζ)

    ϕ = zeros(ComplexF64, grids.r.N, grids.θ.N, ζgrid_size)
    #ϕ = zeros(ComplexF64, size(ϕ_ffs))

    ϕft = fft(ϕ_recon, [2])

    ζgrid = LinRange(0, 2π, ζgrid_size+1)[1:end-1]

    _, _, _, nlist, _= instantiate_grids(grids)

    #this is kind of pointless, given we don't have any way to look
    #at the ζ components yet
    #but future proofing I guess.
    #perhaps we should write our own function for this weird custom ifft.

    for i in 1:grids.r.N
        
        for j in 1:grids.θ.N
            for k in 1:1:length(nlist)

                ϕ[i, j, :] += ϕ_recon[i, j, k] .* exp.(1im * nlist[k] .* ζgrid)
            end
        end
    end

    return ϕ, ϕft

end


function ft_phi(ϕ_recon, grids::FFFGridsT)
    #this is identical regardless of deriv sttaus, not true in general.
    return ϕ_recon, fft(ϕ_recon, [2, 3])

end

function label_mode(ϕft, grids::GridsT)

    _, msize, nsize = phi_size(grids, false)

    rm = zeros(Int64, msize, nsize)
    ϕm = zeros(Float64, msize, nsize)

    for j in 1:msize, k in 1:nsize

        rm[j, k] = argmax(abs.(real.(ϕft[:, j, k])))
        ϕm[j, k] = abs.(real.(ϕft[rm[j, k], j, k]))
    end

    max_mode = argmax(ϕm)

    mlab = mode_label(max_mode[1], grids.θ)
    nlab = mode_label(max_mode[2], grids.ζ)

    rgrid = instantiate_grids(grids)[1] #awful!
    return rgrid[rm[max_mode]], (mlab, nlab)

end

function phi_size(grids::FSSGridsT, ift = false)

    if ift
        θgrid_size = compute_ifft_grid(grids.θ)
        ζgrid_size = compute_ifft_grid(grids.ζ)
    else
        θgrid_size = grids.θ.count
        ζgrid_size = grids.ζ.count
    end

    return grids.r.N, θgrid_size, ζgrid_size
end


function phi_size(grids::FFSGridsT, ift = false)

    if ift
        ζgrid_size = compute_ifft_grid(grids.ζ)
    else
        ζgrid_size = grids.ζ.count
    end

    return grids.r.N, grids.θ.N, ζgrid_size
end


function phi_size(grids::FFFGridsT, ift = false)

    return grids.r.N, grids.θ.N, grids.ζ.N
    
end