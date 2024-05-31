


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