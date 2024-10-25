
#testing a proper q-profiel for damping.

using MID
using Plots; plotlyjs()

#about as big as direct solving can handle
#for this case, fss is much better, 
Nr = 500



geo = GeoParamsT(R0 = 5)
rgrid = rfem_grid(N=Nr)
θgrid = asm_grid(start=1, N=5)
#θgrid = afem_grid(N=10, pf=1)
ζgrid = asm_grid(start=-2, N=1)

grids = init_grids(rgrid, θgrid, ζgrid)

function tae_test_q(r::Float64)
    #should give tae (4, -2) tae at 0.75 with (2, -1) isl at 0.5.
    a = 1.8
    b = 0.8
    q = a + b * r^2
    dq = 2 * b * r
    return q, dq
end


flr = FLRT(δ=-4e-7)
prob = init_problem(q = tae_test_q, geo=geo, flr=flr)


evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=false, target_freq = 0.21, nev=100);

continuum_plot(evals)#, ymax=10)

tae_ind = find_ind(evals, 0.1931409)

display(evals.ω[tae_ind])
println(evals.ω[1:45])
potential_plot(ϕft, grids, tae_ind)