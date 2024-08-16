
#file seems to be ok, hard to tell without better resolution.

function plot_potential(ϕ::Array{ComplexF64}, grids::FSSGridsT, ind=1::Int64; n=nothing, savefile=nothing)

    rgrid, _, mlist, _, _, nlist, _= instantiate_grids(grids)

    #p = plot(r, real.(ϕ[ind, :, 1, n]), label=mlist[1], dpi=600)
    #may want the eigenvalue inside the title?
    p = plot(xlabel=L"r", ylabel=L"\phi", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)

    #this will allow plotting if all eigenfunctions are passed in (i.e. from MID)
    #or if a single is passed in (i.e. from file/MIDParallel)
    if length(size(ϕ)) == 4
        ϕ_plot = ϕ[ind, :, :, :]
    else
        ϕ_plot = ϕ
    end

    if isnothing(n)
        for i in 1:grids.θ.count
            for j in 1:grids.ζ.count
                plot!(rgrid, real.(ϕ_plot[ :, i, j]), label=@sprintf("(%s, %s)", mlist[i], nlist[j]))
            end
        end
    else
        for i in 1:grids.θ.count
            
            plot!(rgrid, real.(ϕ_plot[:, i, n]), label=@sprintf("(%s, %s)", mlist[i], n))
        end
    end

    display(p)

    if !isnothing(savefile)
        savefig(p, savefile)
    end


end


function plot_potential(ϕ::Array{ComplexF64}, grids::FFSGridsT, ind=1::Int64; n=nothing, savefile=nothing)

    #assumes only a single n
    #would be nice if this could do some labelling somehow!
    #but this at least works.

    #mlist = (grids.pmd.start:grids.pmd.incr:grids.pmd.start + grids.pmd.incr * grids.pmd.count)[1:end-1]

    #rgrid = construct_rgrid(grids)

    rgrid, _, _, nlist, _= instantiate_grids(grids)

    #p = plot(r, real.(ϕ[ind, :, 1, n]), label=mlist[1], dpi=600)

    p = plot(xlabel=L"r", ylabel=L"\phi", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)

    if length(size(ϕ)) == 4
        ϕ_plot = ϕ[ind, :, :, :]
    else
        ϕ_plot = ϕ
    end

    #will plot the 1,1 mode twice!
    #this will be way to many modes!
    if isnothing(n)
        for i in 1:grids.θ.N
            #mlab = mod(i-1 + grids.θ.pf, grids.θ.N)
            #think pf is not needed for the labels. v hard to tell though!
            #seems like it is needed for tae modes, but not island modes.
            mlab = mode_label(i, grids.θ)
            for n in 1:grids.ζ.count
                #the order of this is completly cooked, but the labels seem correct.
                
                
                #label should be a function of pf!!!
                #based on simple example with m=1, looks like last ind reps m=0 in some sense, i.e it wraps around?? Not sure how -m's will work for fem2d method.
                plot!(rgrid, real.(ϕ_plot[:, i, n]), label=@sprintf("%s, %s", mlab, nlist[n]))
            end
        end
    else
        for i in 1:grids.θ.N
            #the order of this is completly cooked, but the labels seem correct.
            #mlab = mod(i-1 + grids.θ.pf, grids.θ.N)
            mlab = mode_label(i, grids.θ)
            
            
            #label should be a function of pf!!!
            #based on simple example with m=1, looks like last ind reps m=0 in some sense, i.e it wraps around?? Not sure how -m's will work for fem2d method.
            plot!(rgrid, real.(ϕ_plot[:, i, n]), label=@sprintf("%s, %s", mlab, n))
        end
    end

    display(p)
    if !isnothing(savefile)
        savefig(p, savefile)
    end

    
end




function plot_potential(ϕ::Array{ComplexF64}, grids::FFFGridsT, ind=1::Int64; n=nothing, savefile=nothing)
    #assumes only a single n
    #would be nice if this could do some labelling somehow!
    #but this at least works.

    #mlist = (grids.pmd.start:grids.pmd.incr:grids.pmd.start + grids.pmd.incr * grids.pmd.count)[1:end-1]

    #rgrid = construct_rgrid(grids)

    rgrid, _, _ = instantiate_grids(grids)

    if length(size(ϕ)) == 4
        ϕ_plot = ϕ[ind, :, :, :]
    else
        ϕ_plot = ϕ
    end

    #p = plot(r, real.(ϕ[ind, :, 1, n]), label=mlist[1], dpi=600)

    p = plot(xlabel=L"r", ylabel=L"\phi", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)

    #will plot the 1,1 mode twice!
    #this will be way to many modes!
    if isnothing(n)
        for i in 1:grids.θ.N

            for j in 1:grids.ζ.N
            
                mlab = mode_label(i, grids.θ)
                nlab = mode_label(j, grids.ζ)
                #label should be a function of pf!!!
                #based on simple example with m=1, looks like last ind reps m=0 in some sense, i.e it wraps around?? Not sure how -m's will work for fem2d method.
                plot!(rgrid, real.(ϕ_plot[:, i, j]), label=@sprintf("%s, %s", mlab, nlab))
            end
        end
    else
        for i in 1:grids.θ.N
            
            for j in 1:grids.ζ.N
            
                mlab = mode_label(i, grids.θ)
                nlab = mode_label(j, grids.ζ)
                #i.e only plot desired n.
                if nlab==n
                    #label should be a function of pf!!!
                    #based on simple example with m=1, looks like last ind reps m=0 in some sense, i.e it wraps around?? Not sure how -m's will work for fem2d method.
                    plot!(rgrid, real.(ϕ_plot[:, i, j]), label=@sprintf("%s, %s", mlab, nlab))
                end
            end
        end
    end

    display(p)
    if !isnothing(savefile)
        savefig(p, savefile)
    end

    
end



function find_ind(evals::EvalsT, val)

    return argmin(abs.(evals.ω .-val))
end