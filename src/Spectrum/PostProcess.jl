"""
    function post_process(evals::AbstractArray, efuncs::Array{ComplexF64}, grids::FSSGridsT, geo::GeoParamsT)

Proccesses the eigenvalues and eigenfunctions into useful forms. Returns the an EvalsT struct, containing the normalised eigenvalues and data needed to reconstruct the continuum. Also returns the potential realigned with the 3d grid, with the angular variables not fourier transformed, and with the angular variabels fourier transformed.
"""
function post_process(evals::AbstractArray, efuncs::Array{ComplexF64}, grids::FSSGridsT, geo::GeoParamsT)

    ω = geo.R0 .* sqrt.(evals)

    ϕft = reconstruct_phi(efuncs, length(evals), grids)

    ϕ = zeros(ComplexF64, size(ϕft))

    for i in 1:grids.r.N


        ϕ[:, i, :, :] = ifft(ϕft[:, i, :, :], [2, 3])

    end


    rgrid, _, mlist, _, _, nlist, _= instantiate_grids(grids)


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

    ϕ_ffs = reconstruct_phi(efuncs, length(evals), grids)


    ϕft = zeros(ComplexF64, size(ϕ_ffs))
    ϕ = zeros(ComplexF64, size(ϕ_ffs))

    
    for i in 1:grids.r.N

        for n in 1:grids.ζ.count
        #hopefully the 2d case works as expected
            ϕft[:, i, :, n] = fft(ϕ_ffs[:, i, :, n], [2])
        end

        for j in 1:grids.θ.N
            ϕ[:, i, j, :] = ifft(ϕ_ffs[:, i, j, :], [2])
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