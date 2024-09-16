
using MID
using FFTW
using Plots; plotlyjs()


Nr = 30
Nθ = 6
Nζ = 2
rgrid = init_fem_grid(N=Nr);
θgrid = init_fem_grid(N=Nθ, pf=2);
ζgrid = init_fem_grid(N=Nζ, pf=-2);
#θgrid = init_sm_grid(start=2, count = 2)
#ζgrid = init_sm_grid(start=-2, count=1);

grids = init_grids(rgrid, θgrid, ζgrid);

geo = GeoParamsT(R0=1000.0)

prob = init_problem(q=island_mode_q, geo=geo); 

isl = ContIslandT(m0=3, n0=2, q0=3/2, qp=1.6, A=1e-4, ψ0=0.125)

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true, deriv=true); 


plot_continuum(evals)


ind = find_ind(evals, 0.383)

plot_potential(ϕft[:, :, :, :, 1], grids, ind)


#need some kind of island gridsT tbh.
ϕisl = tor_to_isl(100, 30, 10, ϕ[ind, :, :, :, :, :], grids, isl);

ϕislfft = fft(ϕisl, [2, 3]);

κgrid = LinRange(0, 2, 100)

p = plot()
for i in 1:30

    #plot!(κgrid, real.(ϕislfft[:, i, 1]))
    plot!(κgrid, real.(ϕisl[:, i, 1]))
end
display(p)

#ok so seems to be at least mildly working...
display(ϕislfft[:, 1, 1])