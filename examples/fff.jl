
using MID
using Plots; plotlyjs()

#Integration has been massivly sped up, but this is still slow af, probbaly requires multiproc.
#~90% of the time is spent numerically integrating. Wonder if there is anything we can do???
Nr = 30
Nθ = 6
Nζ = 2
rgrid = rfem_grid(N=Nr, gp=5);
θgrid = afem_grid(N=Nθ, pf=2, gp=5);
ζgrid = afem_grid(N=Nζ, pf=-2, gp=5);

grids = init_grids(rgrid, θgrid, ζgrid);

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo); 

#display(collect(MID.inst_grid(θgrid)))

#W, I = construct(prob, grids);

#display(Matrix(W)[5:8, 5:8])
#display(Matrix(I)[5:8, 5:8])

#println(Matrix(W)[41, :])

#display(size(Matrix(W)))

#with @views. 22.907080 seconds (8.07 M allocations: 721.379 MiB, 1.36% gc time)
#fk load more allocations and gc without views.
#outrageous spead up shifting the ϕ[:, test, :, ...] to ϕ[:, testr, testθ, :, :]

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true); 

continuum_plot(evals)


println((evals.ω  ./ 10) .^2)


ind = find_ind(evals, 0.383)
#ind = 348
println((evals.ω[55:60] ./ 10) .^2)
display((evals.ω[ind] / 10)^2)
potential_plot(ϕft, grids, ind+3)


contour_plot(ϕ, grids, ind=ind)
using FFTW

ft_test = [1.0 + 0.1*1im, 2.0  + 0.2*1im, 3.0 + 0.3im, 4.0 + 0.4im, 5.0 + 0.5im]

fft(ft_test)