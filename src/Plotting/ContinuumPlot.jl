

function plot_continuum(; ω, grids, filename=nothing, n=1)
    #colours of this are cooked, ie gaps dont flip like they should.
    #rotating it gives way to much whitespace between axis label and tickmarks.
    p = scatter(xlabel=L"r", ylabel=L"\frac{\omega  R_0}{v_A}", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)#, guidefontvalign=:hcentre)
    mlist = (grids.pmd.start:grids.pmd.incr:grids.pmd.start + grids.pmd.incr * grids.pmd.count)[1:end-1]
    rgrid = collect(LinRange(0, 1, grids.rd.N))[2:end]
    for (i, m) in enumerate(mlist)
        scatter!(rgrid, ω[:, i], label=@sprintf("(%d, %d)", m, n))
    end
    
    display(p)

    if !isnothing(filename)
        savefig(p, filename)
    end

end




#this probably shouldn't be here!
#note that this method of reconstructing has a few flaws,
#notably that when m=±1, the linear behaviour at r=0 tends to dominate the singularity, so these modes are 
#often places at the wrong radial point. This problem is reduced when the grid resolution goes up as the singularity becomes
#less smooth
#there also seems to be a problem when ω changes sign, in that case we get efuncs with a large width, and they get placed a bit wrong
#this is probably not a huge issue, just annoying that the reconstruction is not perfect
#this is also reduced as the resolution goes up.
#think a better method would be to find the point of greatest change, but that seems more difficult.

function reconstruct_continuum(; ω, ϕ, grids, ymin=-0.05, ymax=1.05, filename=nothing)
    #assumes only 2 modes atm.
    omdata = zeros(length(ω)) 
    rdata = zeros(length(ω))
    #col = zeros(length(ω))

    col = Tuple{Int, Int}[]

    #rm = zeros(Int, pmd.count, tmd.count)
    #ϕm = zeros(pmd.count, tmd.count)
    labels = zeros(grids.pmd.count, grids.tmd.count)

    #may be better to use other modules and call spectral_grid or whatever
    mlist = (grids.pmd.start:grids.pmd.incr:grids.pmd.start + grids.pmd.incr * grids.pmd.count)[1:end-1]
    nlist = (grids.tmd.start:grids.tmd.incr:grids.tmd.start + grids.tmd.incr * grids.tmd.count)[1:end-1]
    
    rgrid = construct_rgrid(grids)
    #labels = [collect(pmd.start:pmd.start+pmd.count), collect(tmd.start:tmd.start+tmd.count)]

    for i in 1:1:length(ω)
        #still assumes only a single n value!
        rm = zeros(Int, grids.pmd.count, grids.tmd.count)
        ϕm = zeros(grids.pmd.count, grids.tmd.count)
        for j in 1:grids.pmd.count

            for k in 1:grids.tmd.count
                rm[j, k] = argmax(abs.(real.(ϕ[i, :, j, k])))
                ϕm[j, k] = abs.(real.(ϕ[i, rm[j, k], j, k]))
            end
        end

        max_mode = argmax(ϕm)

        #display(max_mode)
        #display(max_mode[2])
        #display(rm[max_mode])

        rdata[i] = rgrid[rm[max_mode]]
        #col[i] = (mlist[max_mode[1]], nlist[max_mode[2]])
        push!(col, (mlist[max_mode[1]], nlist[max_mode[2]]))
        omdata[i] = abs.(ω[i]) #already normalised now!
        #omdata[i] = sqrt.(abs.(ω[i])) * R0


    end

    p = scatter(rdata, omdata, group=col, ylimits=(ymin, ymax), xlabel=L"r", ylabel=L"\frac{\omega  R_0}{v_A}", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)
    #return rdata, omdata, col
    display(p)
    if !isnothing(filename)
        savefig(p, filename)
    end
end