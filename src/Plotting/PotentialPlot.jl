"""
    plot_potential(ϕ::Array{ComplexF64}, grids::FSSGridsT, ind=1::Int64; n=nothing, savefile=nothing)


Plots the potential as a function of radius showing the mode structure. Expects the fourier transformed potential.
Can set n to only view poloidal modes of a specific toroidal mode. ind is the index of the eigenfunction, used when an array of eigenfunctions is passed in, and is ignored if only a single eigenfunction is passed in.
"""
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


"""
    plot_potential(ϕ::Array{ComplexF64}, grids::FFSGridsT, ind=1::Int64; n=nothing, savefile=nothing)


Plots the potential as a function of radius showing the mode structure. Expects the fourier transformed potential.
Can set n to only view poloidal modes of a specific toroidal mode. ind is the index of the eigenfunction, used when an array of eigenfunctions is passed in, and is ignored if only a single eigenfunction is passed in.
"""
function plot_potential(ϕ::Array{ComplexF64}, grids::FFSGridsT, ind=1::Int64; n=nothing, savefile=nothing)


    rgrid, _, _, nlist, _= instantiate_grids(grids)


    p = plot(xlabel=L"r", ylabel=L"\phi", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)

    if length(size(ϕ)) == 4
        ϕ_plot = ϕ[ind, :, :, :]
    else
        ϕ_plot = ϕ
    end

    if isnothing(n)
        for i in 1:grids.θ.N

            mlab = mode_label(i, grids.θ)
            for n in 1:grids.ζ.count

                plot!(rgrid, real.(ϕ_plot[:, i, n]), label=@sprintf("%s, %s", mlab, nlist[n]))
            end
        end
    else
        for i in 1:grids.θ.N

            mlab = mode_label(i, grids.θ)
            
        
            plot!(rgrid, real.(ϕ_plot[:, i, n]), label=@sprintf("%s, %s", mlab, n))
        end
    end

    display(p)
    if !isnothing(savefile)
        savefig(p, savefile)
    end

    
end



"""
    plot_potential(ϕ::Array{ComplexF64}, grids::FFFGridsT, ind=1::Int64; n=nothing, savefile=nothing)


Plots the potential as a function of radius showing the mode structure. Expects the fourier transformed potential.
Can set n to only view poloidal modes of a specific toroidal mode. ind is the index of the eigenfunction, used when an array of eigenfunctions is passed in, and is ignored if only a single eigenfunction is passed in.
"""
function plot_potential(ϕ::Array{ComplexF64}, grids::FFFGridsT, ind=1::Int64; n=nothing, savefile=nothing)
   

    rgrid, _, _ = instantiate_grids(grids)

    if length(size(ϕ)) == 4
        ϕ_plot = ϕ[ind, :, :, :]
    else
        ϕ_plot = ϕ
    end


    p = plot(xlabel=L"r", ylabel=L"\phi", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)

    #will plot the 1,1 mode twice!
    #this will be way to many modes!
    if isnothing(n)
        for i in 1:grids.θ.N

            for j in 1:grids.ζ.N
            
                mlab = mode_label(i, grids.θ)
                nlab = mode_label(j, grids.ζ)

                plot!(rgrid, real.(ϕ_plot[:, i, j]), label=@sprintf("%s, %s", mlab, nlab))
            end
        end
    else
        for i in 1:grids.θ.N
            
            for j in 1:grids.ζ.N
            
                mlab = mode_label(i, grids.θ)
                nlab = mode_label(j, grids.ζ)
                if nlab==n

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
