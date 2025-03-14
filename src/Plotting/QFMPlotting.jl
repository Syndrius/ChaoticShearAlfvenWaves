

#this should probably be in plotting eventually!
function plot_surfs(surfs::Array{QFMSurfaceT})
    #expects an array of surface structs.

    display(surfs)

    p = plot(xlabel=L"r", ylabel=L"\phi", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)

    for i in surfs

        s, t = convert_surf(i)

        plot!(t, s)

    end

    display(p)


end