"""
    continuum_plot(evals::EvalsT; savefile=nothing, n=nothing, ymin=-0.05, ymax=1.05)

Plots the conntinuum based on an EvalsT struct, produced from compute_spectrum.
"""
function continuum_plot(evals::EvalsT; savefile=nothing, n=nothing, ymin=-0.05, ymax=1.05)

    if isnothing(n)
        p = scatter(evals.r, real.(evals.ω), group=evals.modelabs, xlabel=L"r", ylabel=L"\frac{\omega  R_0}{v_A}", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10, ylimits=(ymin, ymax))

        display(p)
    else
        #think we assume it is just one value for now.
        #seems like this should be doable in one.
        #bit clunky but does work.
        ns = [x[2] for x in evals.modelabs]
        mls = evals.modelabs[ns .== n]
        rs = evals.r[ns .== n]

        oms = evals.ω[ns .== n]

        p = scatter(rs, real.(oms), group=mls, xlabel=L"r", ylabel=L"\frac{\omega  R_0}{v_A}", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10, ylimits=(ymin, ymax))

        display(p)
    end

    if !isnothing(savefile)
        savefig(p, savefile)
    end
end



"""
    continuum_plot(ω, grids::ContGridsT; savefile=nothing, ymin=-0.05, ymax=1.05)

Plots the continuum produces by compute_continuum().
"""

#this is a shit function.
#don't think this will ever work.
#only works if we do continuum with the perN function
#otherwise this is useless.
function continuum_plot(ω, grids::ContGridsT; savefile=nothing, ymin=-0.05, ymax=1.05)

    #note that this does not label m's, not sure if it is possible
    #could pair this with an analytical version for cylinder??

    p = scatter(ylimits=(ymin, ymax))
    

    rgrid = inst_grids(grids)[1]

    mlist = mode_list(grids.θ)
    nlist = mode_list(grids.ζ)
    if rgrid[1]==0
        #same condition used when recontsructing
        rgrid = rgrid[2:end]
    end
    rgrid_plot = repeat(rgrid, 1, size(ω)[2])
    for (i, n) in enumerate(nlist)
        scatter!(vcat(rgrid_plot...), vcat(ω[:, :, i]...), label=@sprintf("n=%d", n))
    end

    display(p)

    if !isnothing(savefile)
        savefig(p, savefile)
    end

end
