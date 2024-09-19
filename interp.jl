
#just using inbuilt interpolation instead of hermite to see if it is better.


Nr = 30
Nθ = 6
Nζ = 2
rgrid = init_fem_grid(N=Nr);
θgrid = init_fem_grid(N=Nθ, pf=2);
ζgrid = init_fem_grid(N=Nζ, pf=-2);
#θgrid = init_sm_grid(start=2, count = 2)
#ζgrid = init_sm_grid(start=-2, count=1);

grids = init_grids(rgrid, θgrid, ζgrid);

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo); 


#with @views. 22.907080 seconds (8.07 M allocations: 721.379 MiB, 1.36% gc time)
#fk load more allocations and gc without views.
#outrageous spead up shifting the ϕ[:, test, :, ...] to ϕ[:, testr, testθ, :, :]

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true, deriv=false); 


plot_continuum(evals)



ind = find_ind(evals, 0.383)
#ind = 348
plot_potential(ϕft, grids, ind)


contour_plot(ϕ, grids, ind=ind)

test_fn = ϕ[ind, :, :, :];


rg, θg, ζg = instantiate_grids(grids);

itp = interpolate((rg, collect(θg), collect(ζg)), test_fn, Gridded(Linear(Periodic())))

ext = extrapolate(itp, Periodic())

Nrtest = 100
Nθtest = 30
Nζtest = 15

rtest = LinRange(0, 1, Nrtest)
θtest = LinRange(0, 2π, Nθtest+1)[1:end-1]
ζtest = LinRange(0, 2π, Nζtest+1)[1:end-1]

ϕ_int = zeros(ComplexF64, Nrtest, Nθtest, Nζtest);
for (i, r) in enumerate(rtest), (j, θ) in enumerate(θtest), (k, ζ) in enumerate(ζtest)

    ϕ_int[i, j, k] = ext(r, θ, ζ)


end

contourf(θtest, rtest, real.(ϕ_int[:, :, 1]), color=:turbo, levels=100)

ϕ_int_fft = fft(ϕ_int, [2, 3]);

p = plot()
for i in 1:Nθtest
    plot!(rtest, real.(ϕ_int_fft[:, i, 1]))
    #plot!(rtest, real.(ϕ_int[:, i, 1]))
end
display(p)