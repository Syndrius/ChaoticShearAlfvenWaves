


#this is much more important to have the docstrings stuff, always forget these args.
#plots all the poloidal modes of the eigenfunction for a given n.
function plot_potential(; r, ϕ, ind, pmd, n, save=false, savename=nothing)

    #assumes only a single n
    #would be nice if this could do some labelling somehow!
    #but this at least works.

    mlist = (pmd.start:pmd.incr:pmd.start + pmd.incr * pmd.count)[1:end-1]

    p = plot(r, real.(ϕ[ind, :, 1, n]), label=mlist[1], dpi=500)

    #will plot the 1,1 mode twice!
    for i in 2:pmd.count
        
        plot!(r, real.(ϕ[ind, :, i, n]), label=mlist[i])
    end

    if save

        savefig(p, savename)
    end

    display(p)
end


#
function find_ind(ω, val)

    return argmin(abs.(ω.-val))
end