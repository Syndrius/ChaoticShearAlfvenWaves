


#this is much more important to have the docstrings stuff, always forget these args.
#plots all the poloidal modes of the eigenfunction for a given n.
function plot_potential(ϕ, grids::FSSGridsT, ind, n=1, filename=nothing)

    #assumes only a single n
    #would be nice if this could do some labelling somehow!
    #but this at least works.

    rgrid, _, mlist, _, _, _, _= instantiate_grids(grids)

    #p = plot(r, real.(ϕ[ind, :, 1, n]), label=mlist[1], dpi=600)

    p = plot(xlabel=L"r", ylabel=L"\phi", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)

    #will plot the 1,1 mode twice!
    for i in 1:grids.θ.count
        
        plot!(rgrid, real.(ϕ[ind, :, i, n]), label=@sprintf("m=%s", mlist[i]))
    end

    display(p)
    if !isnothing(filename)
        savefig(p, filename)
    end

    
end


#requires that mode structure has been passed in here!
function plot_potential(ϕms, grids::FFSGridsT, ind, n=1, filename=nothing)

    #assumes only a single n
    #would be nice if this could do some labelling somehow!
    #but this at least works.

    #mlist = (grids.pmd.start:grids.pmd.incr:grids.pmd.start + grids.pmd.incr * grids.pmd.count)[1:end-1]

    #rgrid = construct_rgrid(grids)

    rgrid, _, _, _, _= instantiate_grids(grids)

    #p = plot(r, real.(ϕ[ind, :, 1, n]), label=mlist[1], dpi=600)

    p = plot(xlabel=L"r", ylabel=L"\phi", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)

    #will plot the 1,1 mode twice!
    #this will be way to many modes!
    for i in 1:grids.θ.N
        #the order of this is completly cooked, but the labels seem correct.
        mlab = mod(i-1 + grids.θ.pf, grids.θ.N)
        if mlab > grids.θ.N/2
            mlab = mlab - grids.θ.N
        end
        
        #label should be a function of pf!!!
        #based on simple example with m=1, looks like last ind reps m=0 in some sense, i.e it wraps around?? Not sure how -m's will work for fem2d method.
        plot!(rgrid, real.(ϕms[ind, :, i, n]), label=@sprintf("m=%s", mlab))
    end

    display(p)
    if !isnothing(filename)
        savefig(p, filename)
    end

    
end


#requires that mode structure has been passed in here!
function plot_potential(ϕms, grids::FFFGridsT, ind, n=1, filename=nothing)

    #assumes only a single n
    #would be nice if this could do some labelling somehow!
    #but this at least works.

    #mlist = (grids.pmd.start:grids.pmd.incr:grids.pmd.start + grids.pmd.incr * grids.pmd.count)[1:end-1]

    #rgrid = construct_rgrid(grids)

    rgrid, _, _ = instantiate_grids(grids)

    #p = plot(r, real.(ϕ[ind, :, 1, n]), label=mlist[1], dpi=600)

    p = plot(xlabel=L"r", ylabel=L"\phi", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)

    #will plot the 1,1 mode twice!
    #this will be way to many modes!
    for i in 1:grids.θ.N
        
        mlab = mod(i-1 + grids.θ.pf, grids.θ.N)
        if mlab > grids.θ.N/2
            mlab = mlab - grids.θ.N
        end
        #label should be a function of pf!!!
        #based on simple example with m=1, looks like last ind reps m=0 in some sense, i.e it wraps around?? Not sure how -m's will work for fem2d method.
        plot!(rgrid, real.(ϕms[ind, :, i, n]), label=@sprintf("m=%s", mlab))
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


function plot_phi_surface(ϕ, grids::FFSGridsT, ind, n=1, filename=nothing)
    #not sure how to treat n in this case


    rgrid, θgrid, _, _, _ = instantiate_grids(grids)   


    #currently construct only does a specific n.
    p = surface(θgrid, rgrid, real.(ϕ[ind, :, :, n]))
    display(p)
    if !isnothing(filename)
        savefig(p, filename)
    end

end


function plot_phi_surface(ϕ, grids::FFFGridsT, ind, ζ=1, filename=nothing)
    #pass in a zeta index, ideally could map that to a position or something.


    rgrid, θgrid, _ = instantiate_grids(grids)   


    #currently construct only does a specific n.
    p = surface(θgrid, rgrid, real.(ϕ[ind, :, :, ζ]))
    display(p)
    if !isnothing(filename)
        savefig(p, filename)
    end

end


#note requires that phi surface has been called.
function plot_phi_surface(ϕsur, grids::FSSGridsT, ind, filename=nothing)

    rgrid, Nθ, mlist, θgrid, _, _, _ = instantiate_grids(grids)

   
    if Nθ < 50
        Nθ = 50
        θgrid = LinRange(0, 2π, 50)
    end

    #currently construct only does a specific n.
    p = surface(θgrid, rgrid, ϕsur[ind, :, :])
    display(p)
    if !isnothing(filename)
        savefig(p, filename)
    end

end



function construct_surface(ϕ, nevals, grids, ζ)

    #slightly less arbitrary way of doing the θ stuff?
    #this is hopeless when we only use a few modes.
    rgrid, Nθ, mlist, θgrid, _, nlist, _ = instantiate_grids(grids)

    #this is a stupid condition but makes it look smoother so idk
    if Nθ < 50
        Nθ = 50
        θgrid = LinRange(0, 2π, 50)
    end
    #nθ = 50 #this is arbitrary but I don't think it should be!
    #maybe we need to inverse the fourier transform? Maybe that is what we are kind of doing?

    #θgrid = LinRange(0, 2π, nθ) # I fell like this doesn't make sense to make our own arbitrary θgrid
    #seems like it should be aligned with our island or solver or something.

    ϕsur = zeros(nevals, grids.r.N, Nθ) #start with a single n value. can add extra dimension to this for n.

    

    for i in 1:grids.r.N

        for (j, θ) in enumerate(θgrid)

            for (k, m) in enumerate(mlist)

                for (l, n) in enumerate(nlist)

                    #-ve here? perhap?
                    ϕsur[:, i, j] += @. real(ϕ[:, i, k, l] * exp(1im * (m * θ + n * ζ)))
                end
            end
        end
    end

    return ϕsur

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


#probably just for ffs atm.
function contour_plot(ϕ, grids, ind; ymin=nothing, ymax=nothing, filename=nothing)

    rgrid, θgrid, _, _, _ = instantiate_grids(grids)
    z = zeros(Float64, grids.r.N, grids.θ.N)
    for n in grids.ζ.count
        z += ϕ[ind, :, :, n]
    end

    #if nothing, this still work, just gives warning.
    contour(θgrid, rgrid, real.(z), ylimits=(ymin, ymax))

end