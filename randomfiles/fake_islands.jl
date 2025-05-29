
#making fake island modes for testing the mapping

using MID
using MIDViz
using Plots; plotlyjs()
#%%

function inside_island_mode(κ, ᾱ, φ)

    m=2
    if κ > 1
        return 0
    end
    #return exp(-(κ-0.7)^2/0.01) * exp(1im * m * ᾱ)
    return (-4*(κ-0.5)^2+1)*exp(1im * m * ᾱ)
    #return (-4*(κ-0.5)^2+1)*cos(m * ᾱ)
end


function outside_island_mode(κ, ᾱ, φ)
    m=2
    if κ < 1
        return 0
    end
    #assumes κ goes to 8.
    return exp(-(κ-6.7)^2/0.01) * exp(1im * m * ᾱ)
    #return exp(-(κ-6.7)^2/0.01) * cos(m * ᾱ)
end

function create_island_mode(κgrid, ᾱgrid, φgrid, inside)

    ϕ_isl = zeros(ComplexF64, length(κgrid), length(ᾱgrid), length(φgrid))

    if inside

        for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, φ) in enumerate(φgrid)
            ϕ_isl[i, j, k] = inside_island_mode(κ, ᾱ, φ)
        end
    else
        #κgrid = LinRange(0, 5, Nκ)
        for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, φ) in enumerate(φgrid)
            ϕ_isl[i, j, k] = outside_island_mode(κ, ᾱ, φ)
        end

    end
    return ϕ_isl

end

function map_plot(ϕ, mapgrids)
    x1grid, _, _ = MIDIslands.Mapping.inst_grids(mapgrids)

    p = plot()
    for i in 1:size(ϕ)[2]
        plot!(x1grid, real.(ϕ[:, i, 1]))
    end
    display(p)
end
function map_contour(ϕ, mapgrids)
    x1grid, x2grid, _ = MIDIslands.Mapping.inst_grids(mapgrids)
    #tis a problem that this needs to occur.
    #says that something is wrong with our transformation.
    #this is only needed for the final case, otherwise it cooked it lol
    #coordinate transformation needs to be perfected.
    contourf(x2grid, sqrt.(x1grid), real.(ϕ[:, :, 1]))
end
#%%

isl = init_island(m0=2, n0=-1, w=0.03, r0=0.5, qp=2.0)
isl = MID.Geometry.inst_island(isl)
#%%
Nκ = 500
Nᾱ = 100
Nφ = 10
κf = LinRange(0, 0.999, Nκ)
ᾱf = LinRange(0, 2π, Nᾱ+1)[1:end-1]
φf = LinRange(0, 2π, Nφ+1)[1:end-1]

#idealised island mode, created in island geometry.
ϕf = create_island_mode(κf, ᾱf, φf, true);
contourf(ᾱf, κf, real.(ϕf[:, :, 1]))
#%%

#now we transform into toroidal coordinates
#this being a mapped grid is going to cause problemos
κgrid = init_grid(type=:rf, N=Nκ, stop=0.999)
ᾱgrid = init_grid(type=:af, N=Nᾱ) #probs no pf for now, shouldn't matter
φgrid = init_grid(type=:af, N=Nφ) #probs no pf for now, shouldn't matter
isl_grids = init_grids(κgrid, ᾱgrid, φgrid)
Nr = 500
Nθ = 100
Nζ = 10
rgrid = init_grid(type=:rf, N=Nr)
θgrid = init_grid(type=:af, N=Nθ)
ζgrid = init_grid(type=:af, N=Nζ)
tor_grids = init_grids(rgrid, θgrid, ζgrid)
#tor_mgrids = MIDIslands.Mapping.MapGridsT(Nx1=Nr, x1max=1.0, Nx2=Nθ, Nx3=Nζ)


ϕtor = map_isl_to_tor(tor_grids, ϕf, isl_grids, isl);
#%%

#seems to be ok, not convinced tbh
potential_plot(ϕtor, tor_grids)
MIDViz.Plotting.contour_plot(ϕtor, tor_grids)
#%%
#isl_mgrids = MIDIslands.Mapping.MapGridsT(Nx1=500, x1max=1.0, Nx2=40, Nx3=5)

ϕisl_in, ϕisl_outp, ϕisl_outm = map_tor_to_isl(isl_grids, ϕtor, tor_grids, isl);
#%%
potential_plot(ϕisl_in, isl_grids)
MIDViz.Plotting.contour_plot(ϕisl_in, isl_grids)

#%%
#testing the periodicity of interpolations
using Interpolations

N = 45
x = collect(LinRange(0, 2π, N+1)[1:end-1])
y = collect(LinRange(0, 2π, N+1)[1:end-1])
f = cos.(x) .+ sin.(y') 


S = interpolate((x, y), f, (Gridded(Linear(Periodic(OnCell()))), Gridded(Linear(Periodic(OnCell())))))
Sext = extrapolate(S, Periodic());
#%%
N2 = 100
xi = LinRange(0, 2π, N2+1)[1:end-1]
yi = LinRange(0, 2π, N2+1)[1:end-1]
fi = zeros(N2, N2);

#S = interpolate((xi, yi), f, (Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))))

for (i, xval) in enumerate(xi), (j, yval) in enumerate(yi)
    fi[i, j] = Sext(xval, yval)
end

#%%

#this package is fkn awful
S = interpolate( cos.(x), BSpline(Cubic(Periodic())));
Sext = extrapolate(S, Periodic())

plot(xi, Sext(xi))
plot!(xi, cos.(xi))


#plot(xi, fi[:, 1], legend=false)
