"""
    contour_plot(ϕ, grids::FFSGridsT, ind, ζ=1; savefile=nothing)

Plots the contours for a given ζ slice. Expects the un-fourier transformed solution.
"""
function contour_plot(ϕ, grids::FFSGridsT, ζ=1; ind, savefile=nothing)

    if length(size(ϕ)) == 4
        ϕ_plot = ϕ[ind, :, :, :]
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
function contour_plot(ϕ, grids::FFFGridsT, ζ=1; ind, savefile=nothing)

    if length(size(ϕ)) == 4
        ϕ_plot = ϕ[ind, :, :, :]
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
function contour_plot(ϕ, grids::FSSGridsT, ζ=1; ind, savefile=nothing)

    if length(size(ϕ)) == 4
        ϕ_plot = ϕ[ind, :, :, :]
    else
        ϕ_plot = ϕ
    end
    rgrid, _, _, _, _, _, _ = instantiate_grids(grids)

    #adds back the periodicity.
    z = zeros(ComplexF64, grids.r.N, grids.θ.count+1)

    z[:, 1:end-1] = ϕ_plot[:, :, ζ]
    z[:, end] = ϕ_plot[ :, 1, ζ]

    θgrid = range(0, 2π, grids.θ.count+1)
    p = contourf(θgrid, rgrid, real.(z), levels=100, color=:turbo)

    display(p)

    if !isnothing(savefile)
        savefig(p, savefile)
    end

end