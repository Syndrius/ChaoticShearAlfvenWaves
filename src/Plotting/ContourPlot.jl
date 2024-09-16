"""
    contour_plot(ϕ, grids::FFSGridsT, ind, ζ=1; savefile=nothing)

Plots the contours for a given ζ slice. Expects the un-fourier transformed solution.
"""
function contour_plot(ϕ, grids::FFSGridsT, ζ=1; ind=1, savefile=nothing)

    #there are more cases now big rip
    if length(size(ϕ)) == 5 #all solutions and derivs case
        ϕ_plot = ϕ[ind, :, :, :, 1]
    elseif length(size(ϕ)) == 4
        #probably the derivative case, with both conditions this should be fine right???
        if size(ϕ)[end] == 8 && size(ϕ)[1] == grids.r.N 
            ϕ_plot = ϕ[:, :, :, 1]
        else
            ϕ_plot = ϕ[ind, :, :, :]
        end

    else
        ϕ_plot = ϕ
    end
    rgrid, θgrid, _, _, _ = instantiate_grids(grids)

    #adds back the periodicity.
    z = zeros(ComplexF64, grids.r.N, grids.θ.N+1)

    z[:, 1:end-1] = ϕ_plot[ :, :, ζ]
    z[:, end] = ϕ_plot[ :, 1, ζ]

    θgrid = range(0, 2π, grids.θ.N+1)
    p = contourf(θgrid, rgrid, real.(z), levels=100, color=:turbo)

    display(p)

    if !isnothing(savefile)
        savefig(p, savefile)
    end

end


"""
    contour_plot(ϕ, grids::FFFGridsT, ind, ζ=1; savefile=nothing)

Plots the contours for a given ζ slice. Expects the un-fourier transformed solution.
"""
function contour_plot(ϕ, grids::FFFGridsT, ζ=1; ind=1, savefile=nothing)

    #there are more cases now big rip
    if length(size(ϕ)) == 5 #all solutions and derivs case
        ϕ_plot = ϕ[ind, :, :, :, 1]
    elseif length(size(ϕ)) == 4
        #probably the derivative case, with both conditions this should be fine right???
        if size(ϕ)[end] == 8 && size(ϕ)[1] == grids.r.N 
            ϕ_plot = ϕ[:, :, :, 1]
        else
            ϕ_plot = ϕ[ind, :, :, :]
        end

    else
        ϕ_plot = ϕ
    end
    rgrid, θgrid, _ = instantiate_grids(grids)

    #adds back the periodicity.
    z = zeros(ComplexF64, grids.r.N, grids.θ.N+1)

    z[:, 1:end-1] = ϕ_plot[ :, :, ζ]
    z[:, end] = ϕ_plot[ :, 1, ζ]

    θgrid = range(0, 2π, grids.θ.N+1)

    p = contourf(θgrid, rgrid, real.(z), levels=100, color=:turbo)

    display(p)

    if !isnothing(savefile)
        savefig(p, savefile)
    end

end


"""
    contour_plot(ϕ, grids::FSSGridsT, ind, ζ=1; savefile=nothing)

Plots the contours for a given ζ slice. Expects the un-fourier transformed solution.
"""
function contour_plot(ϕ, grids::FSSGridsT, ζ=1; ind=1, savefile=nothing)

    #there are more cases now big rip
    if length(size(ϕ)) == 5 #all solutions and derivs case
        ϕ_plot = ϕ[ind, :, :, :, 1]
    elseif length(size(ϕ)) == 4
        #probably the derivative case, with both conditions this should be fine right???
        if size(ϕ)[end] == 8 && size(ϕ)[1] == grids.r.N 
            ϕ_plot = ϕ[:, :, :, 1]
        else
            ϕ_plot = ϕ[ind, :, :, :]
        end

    else
        ϕ_plot = ϕ
    end
    rgrid, _, _, _, _, _, _ = instantiate_grids(grids)

    θgrid_size = compute_ifft_grid(grids.θ)
    #adds back the periodicity.
    z = zeros(ComplexF64, grids.r.N, θgrid_size+1)

    z[:, 1:end-1] = ϕ_plot[:, :, ζ]
    z[:, end] = ϕ_plot[ :, 1, ζ]

    #θgrid = range(0, 2π, grids.θ.count+1)
    θgrid = range(0, 2π, θgrid_size+1)

    p = contourf(θgrid, rgrid, real.(z), levels=100, color=:turbo)

    display(p)

    if !isnothing(savefile)
        savefig(p, savefile)
    end

end