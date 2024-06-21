


#this is much more important to have the docstrings stuff, always forget these args.
#plots all the poloidal modes of the eigenfunction for a given n.
function plot_potential(; grids, ϕ, ind, n, filename=nothing)

    #assumes only a single n
    #would be nice if this could do some labelling somehow!
    #but this at least works.

    mlist = (grids.pmd.start:grids.pmd.incr:grids.pmd.start + grids.pmd.incr * grids.pmd.count)[1:end-1]

    rgrid = construct_rgrid(grids)

    #p = plot(r, real.(ϕ[ind, :, 1, n]), label=mlist[1], dpi=600)

    p = plot(xlabel=L"r", ylabel=L"\phi", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)

    #will plot the 1,1 mode twice!
    for i in 1:grids.pmd.count
        
        plot!(rgrid, real.(ϕ[ind, :, i, n]), label=@sprintf("m=%s", mlist[i]))
    end

    display(p)
    if !isnothing(filename)
        savefig(p, filename)
    end

    
end


#
function find_ind(ω, val)

    return argmin(abs.(ω.-val))
end


function plot_surface(z, grids, ind)

    rgrid = construct_rgrid(grids)
    nθ, mlist, θgrid = spectral_grid(grids.pmd)

    if nθ < 50
        nθ = 50
        θgrid = LinRange(0, 2π, 50)
    end

    #currently construct only does a specific n.
    surface(θgrid, rgrid, z[ind, :, :])


end



function construct_surface(ϕ, nevals, grids, ζ)

    #slightly less arbitrary way of doing the θ stuff?
    #this is hopeless when we only use a few modes.
    nθ, mlist, θgrid = spectral_grid(grids.pmd)
    _, nlist, _ = spectral_grid(grids.tmd)

    if nθ < 50
        nθ = 50
        θgrid = LinRange(0, 2π, 50)
    end
    #nθ = 50 #this is arbitrary but I don't think it should be!
    #maybe we need to inverse the fourier transform? Maybe that is what we are kind of doing?

    #θgrid = LinRange(0, 2π, nθ) # I fell like this doesn't make sense to make our own arbitrary θgrid
    #seems like it should be aligned with our island or solver or something.

    z = zeros(nevals, grids.rd.N, nθ) #start with a single n value. can add extra dimension to this for n.

    

    for i in 1:grids.rd.N

        for (j, θ) in enumerate(θgrid)

            for (k, m) in enumerate(mlist)

                for (l, n) in enumerate(nlist)

                    z[:, i, j] += @. real(ϕ[:, i, k, l] * exp(1im * (m * θ + n * ζ)))
                end
            end
        end
    end

    return z

end


#plots the potential but sums different n values.
function plot_sum_potential(; grids, ϕ, ind, filename=nothing)

    #assumes only a single n
    #would be nice if this could do some labelling somehow!
    #but this at least works.

    mlist = (grids.pmd.start:grids.pmd.incr:grids.pmd.start + grids.pmd.incr * grids.pmd.count)[1:end-1]

    rgrid = construct_rgrid(grids)

    #p = plot(r, real.(ϕ[ind, :, 1, n]), label=mlist[1], dpi=600)

    p = plot(xlabel=L"r", ylabel=L"\phi", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)

    #will plot the 1,1 mode twice!
    for i in 1:grids.pmd.count
        res = zeros(size(ϕ[ind, :, i, 1]))

        for j in 1:grids.tmd.count
            #not sure what if this is correct, may need fourier exponent of each n part??
            res += real.(ϕ[ind, :, i, j])
        end
        plot!(rgrid, res, label=@sprintf("m=%s", mlist[i]))
        #plot!(rgrid, real.(ϕ[ind, :, i, n]), label=@sprintf("m=%s", mlist[i]))
    end

    display(p)
    if !isnothing(filename)
        savefig(p, filename)
    end

    
end