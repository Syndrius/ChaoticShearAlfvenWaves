
#this file was creating to swap r to psi, which we have decided not to do.

using MID
using Plots; plotlyjs()


isl = MID.Geometry.init_island(w=0.05, m0=3, n0=-2)

isl = MID.Geometry.inst_island(isl, island_mode_q)

display(isl)
#start very small, matrix scales much more extremly
Nr = 50;

geo = GeoParamsT(R0=10.0)

#prob = init_problem(q=MID.MagneticField.generic_island_q, geo=geo, met=MID.Geometry.cylindrical_metric!)#, isl=isl)#, met=no_delta_metric!); 
prob = init_problem(q=Axel_q, geo=geo)#, met=MID.Geometry.cylindrical_metric!)#, isl=isl)

rgrid = rfem_grid(N=Nr)#, start=0.3, stop=1)
#θgrid = asm_grid(N=2, start=2)
θgrid = afem_grid(N=5, pf=2)
ζgrid = asm_grid(N=1, start=-2)
#ζgrid = afem_grid(N=2, pf=-2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)


#@profview evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);

#scatter(cr.r, real.(cr.ω), ylimits=(-0.05, 1.05))
continuum_plot(evals)

#tae_ind = find_ind(evals, 0.3764079)
tae_ind = find_ind(evals, 0.384)

contour_plot(ϕ, grids, ind=tae_ind)

#display(prob)

potential_plot(ϕft, grids, tae_ind)
potential_plot(ϕ, grids, tae_ind)




display(length(methods(MID.MagneticField.generic_island_q)[1].sig.parameters))
display(length(methods(MID.MagneticField.Axel_q)[1].sig.parameters))


function test_sym(a::Symbol)

    if a==:tub
        display("big boi")
    end
end

test_sym(:hello)

using FFTW

x = rand(100)

p = plan_fft(x)

typeof(p) <: FFTW.FFTWPlan

typeof(p)

