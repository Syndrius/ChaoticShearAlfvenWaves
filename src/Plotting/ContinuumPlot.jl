

#only work for FSS atm. not sure if it will work in general tbh!
#this function is cooked and needs to be fixed.
function plot_continuum(; ω, grids::FSSGridsT, filename=nothing, n=1, ymin=-0.05, ymax=1.05)
    #colours of this are cooked, ie gaps dont flip like they should.
    #rotating it gives way to much whitespace between axis label and tickmarks.
    p = scatter(xlabel=L"r", ylabel=L"\frac{\omega  R_0}{v_A}", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10, ylimits=(ymin, ymax))#, guidefontvalign=:hcentre)
    rgrid, _, mlist, _, _, _, _ = instantiate_grids(grids)

    if rgrid[1]==0
        #same condition used when recontsructing
        rgrid = rgrid[2:end]
    end
    #mlist = (grids.pmd.start:grids.pmd.incr:grids.pmd.start + grids.pmd.incr * grids.pmd.count)[1:end-1]
    #rgrid = collect(LinRange(0, 1, grids.rd.N))[2:end]
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

function reconstruct_continuum(ω, ϕ, grids::FSSGridsT, ymin=-0.05, ymax=1.05, filename=nothing)
    #assumes only 2 modes atm.
    omdata = zeros(length(ω)) 
    rdata = zeros(length(ω))
    #col = zeros(length(ω))

    col = Tuple{Int, Int}[]

    #rm = zeros(Int, pmd.count, tmd.count)
    #ϕm = zeros(pmd.count, tmd.count)
    #labels = zeros(grids.pmd.count, grids.tmd.count)

    #may be better to use other modules and call spectral_grid or whatever
    #mlist = (grids.θ.start:grids.θ.incr:grids.θ.start + grids.θ.incr * grids.θ.count)[1:end-1]
    #nlist = (grids.ζ.start:grids.ζ.incr:grids.ζ.start + grids.ζ.incr * grids.ζ.count)[1:end-1]
    
    #rgrid = construct_rgrid(grids)
    rgrid, _, mlist, _, _, nlist, _= instantiate_grids(grids)
    #labels = [collect(pmd.start:pmd.start+pmd.count), collect(tmd.start:tmd.start+tmd.count)]

    for i in 1:1:length(ω)
        #still assumes only a single n value!
        rm = zeros(Int64, grids.θ.count, grids.ζ.count)
        ϕm = zeros(grids.θ.count, grids.ζ.count)
        for j in 1:grids.θ.count

            for k in 1:grids.ζ.count
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



#note this requires that ϕ has been transformed into the mode structure version.
#either here or in the mode structure it might be nice to remove the v large modes perhaps.
function reconstruct_continuum(ω, ϕms, grids::FFSGridsT, filename=nothing, ymin=-0.05, ymax=1.05)
    omdata = zeros(length(ω)) 
    rdata = zeros(length(ω))
    #col = zeros(length(ω))

    #nlist = (grids.tmd.start:grids.tmd.incr:grids.tmd.start + grids.tmd.incr * grids.tmd.count)[1:end-1]

    rgrid, _, _, nlist, _= instantiate_grids(grids)

    col = Tuple{Int, Int}[]

    #wont generalise to clustred grids obvs.
    #rgrid = LinRange(0, 1, grids.rd.N)
    #θgrid = LinRange(0, 2π, grids.θd.N)

    for i in 1:1:length(ω)

        rm = zeros(Int64, grids.θ.N, grids.ζ.count)
        ϕm = zeros(Float64, grids.θ.N, grids.ζ.count)

        for j in 1:grids.θ.N, k in 1:grids.ζ.count

            rm[j, k] = argmax(abs.(real.(ϕms[i, :, j, k])))
            ϕm[j, k] = abs.(real.(ϕms[i, rm[j, k], j, k]))
        end

        max_mode = argmax(ϕm)

        rdata[i] = rgrid[rm[max_mode]]
        #col[i] = (mlist[max_mode[1]], nlist[max_mode[2]])
        push!(col, (max_mode[1], nlist[max_mode[2]]))
        #probably don't need this anymore tbh.
        omdata[i] = abs.(ω[i]) #already normalised now!
    end

    p = scatter(rdata, omdata, group=col, ylimits=(ymin, ymax))#, xlabel=L"r", ylabel=L"\frac{\omega  R_0}{v_A}", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)
    #return rdata, omdata, col
    display(p)
    if !isnothing(filename)
        savefig(p, filename)
    end
end


#note this requires that ϕ has been transformed into the mode structure version.
#either here or in the mode structure it might be nice to remove the v large modes perhaps.
function reconstruct_continuum(ω, ϕms, grids::FFFGridsT, filename=nothing, ymin=-0.05, ymax=1.05)
    omdata = zeros(length(ω)) 
    rdata = zeros(length(ω))
    #col = zeros(length(ω))

    #nlist = (grids.tmd.start:grids.tmd.incr:grids.tmd.start + grids.tmd.incr * grids.tmd.count)[1:end-1]

    rgrid, _, _, = instantiate_grids(grids)

    col = Tuple{Int, Int}[]

    #wont generalise to clustred grids obvs.
    #rgrid = LinRange(0, 1, grids.rd.N)
    #θgrid = LinRange(0, 2π, grids.θd.N)

    for i in 1:1:length(ω)

        rm = zeros(Int64, grids.θ.N, grids.ζ.N)
        ϕm = zeros(Float64, grids.θ.N, grids.ζ.N)

        for j in 1:grids.θ.N, k in 1:grids.ζ.N

            rm[j, k] = argmax(abs.(real.(ϕms[i, :, j, k])))
            ϕm[j, k] = abs.(real.(ϕms[i, rm[j, k], j, k]))
        end

        max_mode = argmax(ϕm)

        rdata[i] = rgrid[rm[max_mode]]
        #col[i] = (mlist[max_mode[1]], nlist[max_mode[2]])
        #not sure how this will handle the different n's tbh!
        #we probably need a bigger version to test this out bh.
        push!(col, (max_mode[1], 1))
        #probably don't need this anymore tbh.
        omdata[i] = abs.(ω[i]) #already normalised now!
    end

    p = scatter(rdata, omdata, group=col, ylimits=(ymin, ymax))#, xlabel=L"r", ylabel=L"\frac{\omega  R_0}{v_A}", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)
    #return rdata, omdata, col
    display(p)
    if !isnothing(filename)
        savefig(p, filename)
    end
end




#this actually seems to work pretty darn well!
#unsure why the two modes have a different sign in this case, seems a wee bit odd. -> with higher res they have the same sign now, seems odd it ever swapped though!
#this seems like it will be quite flawed for arpack solve, as then there will be less of a fourier transform perhaps?
#not sure, but the m=2 modes tructure seems significantly smaller.
#I don't think this should be here tbh!
function mode_structure(ϕ, grids::FFSGridsT)

    ϕms = zeros(ComplexF64, size(ϕ))

    #i guess we fourier transform at each radial point? Perhap?

    #rgrid = LinRange(0, 1, grids.rd.N) #this obvs won't work in general.

    #may be necesary that we only keep the nth largest modes or something tbh.
    #this will get v large as θ gets v large.
    for i in 1:grids.r.N
        #assume single n for now!
        for n in 1:grids.ζ.count
            ϕms[:, i, :, n] = fft(ϕ[:, i, :, n], [2])
        end

    end

    return ϕms
end


function mode_structure(ϕ, grids::FFFGridsT)

    ϕms = zeros(ComplexF64, size(ϕ))

    #i guess we fourier transform at each radial point? Perhap?

    #rgrid = LinRange(0, 1, grids.rd.N) #this obvs won't work in general.

    #may be necesary that we only keep the nth largest modes or something tbh.
    #this will get v large as θ gets v large.
    for i in 1:grids.r.N
        #hopefully the 2d case works as expected
        ϕms[:, i, :, :] = fft(ϕ[:, i, :, :], [2, 3])

    end

    return ϕms
end