

#this should probably be in plotting eventually!
function plot_surfs(surfs::Array{QFMSurfaceT}; filename=nothing)
    #expects an array of surface structs.

    #display(surfs)

    p = plot(xlabel=L"\theta", ylabel=L"r", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)

    for i in surfs

        s, t = convert_surf(i)

        #plot!(mod.(t, 2π), s)
        #scatter!(mod.(t, 2π), s)
        #neither of hte above do what we want!
        #need to fix this eventually!phiphi
        plot!(t, s)

    end


    if !isnothing(filename)
        savefig(p, filename)
    end

    display(p)

end
