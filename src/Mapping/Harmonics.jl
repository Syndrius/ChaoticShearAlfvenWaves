
#this is perhaps an odd place, 
#but this provides extra information based on the harmonic structure of the spectrum

#terrible name but no idea what to call it.
#if this is vaguely succsesful we probably want to move this into Structures/output.jl
struct HarmonicsT
    fwhm :: Array{Float64} #full width half maximum of dominant harmonic
    dd :: Array{Float64} # second derivative at the peak location of the dominant harmonic
    harm2 :: Array{Float64} # relative peak amplitude of the second largest harmonic
    modelabs :: Array{Tuple{Int64, Int64}} # same as evalsT, so we can plot nicely
    x1 :: Array{Float64} # same as evalsT, the radial location of the dominant harmonic for plotting.
    function HarmonicsT(N::Int64)
        new(zeros(N), zeros(N), zeros(N), Array{Tuple{Int64, Int64}}(undef, N), zeros(N))
    end
end


#finds where a 1d function crosses 0.5
#this includes normalising.
function crossings(f, N)

    am = argmax(abs.(real.(f)))
    nf = real.(f) / real(f[am]) .- 0.5

    nflip = 0
    locs = []
    for i in 1:N-1 #ignore start and finish as they are zero. -> now kept for real edge case where spike is at the edge or start
        if sign(real(nf[i])) ≠ sign(real(nf[i+1]))
            nflip += 1
            #returns the index of the point to the left of the crossing
            push!(locs, i)
        end
    end
    return nflip, locs 
end


function harmonic_info(dir::String)

    prob, grids, _ = inputs_from_file(dir=dir)

    evals = evals_from_file(dir=dir)

    un_inds = load_object(joinpath(dir, "unique_inds.jld2"))

    nevals = length(un_inds)

    ht = HarmonicsT(nevals)

    x1grid = inst_grid(grids.x1)

    x1m = zeros(Int64, grids.x2.N, grids.x3.N)
    intsize = 10
    intg1 = zeros(intsize)
    intg2 = zeros(intsize)
    ϕ1 = zeros(ComplexF64, intsize)
    ϕ2 = zeros(ComplexF64, intsize)

    ϕ, ϕft = PostProcessing.allocate_phi_arrays(grids, deriv=true)
    plan = PostProcessing.create_ft_plan(ϕft, grids)

    ϕm = zeros(grids.x2.N, grids.x3.N)
    ht.modelabs .= evals.modelabs

    for i in 1:nevals

        efunc_read = @sprintf("efunc%05d.hdf5", un_inds[i])
        #unfort doesn't handle complex numbers v well
        efunc_split = load_object(dir*"/efuncs_raw/"*efunc_read)

        #ideally this would be preallocated in some way
        efunc = efunc_split[1, :] .+ efunc_split[2, :] * 1im

        PostProcessing.reconstruct_phi!(efunc, grids, ϕ, ϕft, plan)

        for j in 1:grids.x2.N, k in 1:grids.x3.N
            #note we do not need the derivs for this case~
            x1m[j, k] = argmax(abs.(real.(ϕft[:, j, k, 1])))
            ϕm[j, k] = abs.(real.(ϕft[x1m[j, k], j, k, 1]))
        end
        max_mode = argmax(ϕm)
        max_val = ϕm[max_mode]
        ϕm[max_mode] = 0.0 #so we can find the second
        max2 = argmax(ϕm)

        ht.harm2[i] = ϕm[max2] / max_val
        ht.x1[i] = x1grid[x1m[max_mode]]

        #computes the second derivative at the peak.
        #think this will be completly useless tbh!
        ht.dd[i] = real(hermite_interpolation(ht.x1[i], ϕft[:, max_mode, :], x1grid, 2))

        nflip, locs = crossings(ϕft[:, max_mode, 1], grids.x1.N)
        if length(locs) > 2
            locs = [locs[1], locs[2]]
        end

        intg1 .= LinRange(x1grid[locs[1]], x1grid[locs[1]+1], intsize)
        intg2 .= LinRange(x1grid[locs[2]], x1grid[locs[2]+1], intsize)

        for i in 1:intsize
            ϕ1[i] = hermite_interpolation(intg1[i], ϕft[:, max_mode, :], x1grid)
            ϕ2[i] = hermite_interpolation(intg2[i], ϕft[:, max_mode, :], x1grid)
        end

        ind1 = find_ind(ϕ1, 0.5)
        ind2 = find_ind(ϕ2, 0.5)
        ht.fwhm[i] = intg2[ind2] - intg1[ind1]
    end

    save_object(joinpath(dir, "ht.jld2"), ht)
end

