

#probably just for ffs atm.
function contour_plot(ϕ, grids::FFSGridsT, ind; ymin=nothing, ymax=nothing, filename=nothing)

    rgrid, θgrid, _, _, _ = instantiate_grids(grids)
    θgrid = range(0, 2π, grids.θ.N)
    z = zeros(Float64, grids.r.N, grids.θ.N)
    for n in 1:grids.ζ.count
        z += ϕ[ind, :, :, n]
    end

    #if nothing, this still work, just gives warning.
    contourf(θgrid, rgrid, real.(z), levels=100, color=:turbo,)# ylimits=(ymin, ymax))

end



function contour_plot(ϕ, grids::FSSGridsT, ind; ymin=nothing, ymax=nothing, filename=nothing)
    Nθ = 50
    rgrid, _, mlist, _, _, _, _ = instantiate_grids(grids)
    θgrid = range(0, 2π, Nθ + 1)[1:end-1]
    z = zeros(Float64, grids.r.N, Nθ)
    for (j, m) in enumerate(mlist)
        for n in 1:grids.ζ.count
            for i in 1:Nθ
                z[:, i] += real(ϕ[ind, :, j, n] .* exp(1im*m *θgrid[i]))
            end
        end
    end

    #if nothing, this still work, just gives warning.
    contour(θgrid, rgrid, real.(z), levels = 50)

end


#this is a stupid function that doesn't belong here, should probbaly be in MIDViz.
function plot_contour_poincare(ϕ, grids::FFSGridsT, ind, rp, θp; ymin=0.0, ymax=1.0, filename=nothing)

    rgrid, θgrid, _, _, _ = instantiate_grids(grids)
    #hack to fix plotting, not sure what we are supposed to do tbh
    θgrid = range(0, 2π, grids.θ.N)
    z = zeros(Float64, grids.r.N, grids.θ.N)
    for n in 1:grids.ζ.count
        z += ϕ[ind, :, :, n]
    end

    ms = 1.3
    Ntraj = size(rp)[1]
    p = scatter(markersize=ms, dpi=600)
    for i in 1:Ntraj
        scatter!(rp[i, :], θp[i, :], markersize=ms, alpha=0.2, label=false)
    end
   

    #if nothing, this still work, just gives warning.
    contour!(θgrid, rgrid, real.(z), levels=50, fill=true, color=:turbo, ylimits=(ymin, ymax))

    display(p)
    if !isnothing(filename)
        savefig(p, filename)
    end

end

#maybe rp for poincare not d?
function plot_contour_poincare(ϕ, grids::FSSGridsT, ind, rd, θd; ymin=nothing, ymax=nothing, filename=nothing)
    Nθ = 50
    rgrid, _, mlist, _, _, _, _ = instantiate_grids(grids)
    θgrid = range(0, 2π, Nθ + 1)[1:end-1]
    z = zeros(Float64, grids.r.N, Nθ)
    for (j, m) in enumerate(mlist)
        for n in 1:grids.ζ.count
            for i in 1:Nθ
                z[:, i] += real(ϕ[ind, :, j, n] .* exp(1im*m *θgrid[i]))
            end
        end
    end

    ms = 1.3
    Ntraj = size(rd)[1]
    p = scatter(markersize=ms, legend=false, dpi=600)
    for i in 1:Ntraj
        scatter!(rd[i, :], θd[i, :], markersize=ms, alpha=0.2)
    end

    #if nothing, this still work, just gives warning.
    contour!(θgrid, rgrid, real.(z), levels = 50, fill=true)

    display(p)
    if !isnothing(filename)
        savefig(p, filename)
    end

end
