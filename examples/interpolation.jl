
using MID
using MIDViz
using Plots
using FFTW
#%%
#Integration has been massivly sped up, but this is still slow af, probbaly requires multiproc.
#~90% of the time is spent numerically integrating. Wonder if there is anything we can do???
Nr = 10
Nθ = 3
Nζ = 2
rgrid = init_grid(type=:rf, N=Nr);
θgrid = init_grid(type=:af, N=Nθ, pf=1);
ζgrid = init_grid(type=:af, N=Nζ, pf=-1);
grids = init_grids(rgrid, θgrid, ζgrid);
#%%

#first define the problem
geo = init_geo(R0=4.0)

prob = init_problem(q=fu_dam_q, geo=geo); 
#%%
solver = init_solver(prob=prob, full_spectrum=true)
#%%

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver, deriv=true);
#%%
#scatter(cr.r, real.(cr.ω), ylimits=(-0.05, 1.05))

continuum_plot(evals, n=-1)


ind = find_ind(evals, 0.3)
#ind = 348
potential_plot(ϕft, grids, ind)
#%%
#this should actually be a function inside MID tbh!
Nmr = 50
Nmθ = 20
Nmζ = 10

mrgrid = init_grid(type=:rf, N=Nmr);
mθgrid = init_grid(type=:af, N=Nmθ);
mζgrid = init_grid(type=:af, N=Nmζ);
mgrids = init_grids(mrgrid, mθgrid, mζgrid);

mrg, mθg, mζg = MID.inst_grids(mgrids)

int_phi = zeros(ComplexF64, Nmr, Nmθ, Nmζ);
rg, θg, ζg = MID.inst_grids(grids)
for k in 1:Nmζ, j in 1:Nmθ, i in 1:Nmr
    int_phi[i, j, k] = MID.Mapping.hermite_interpolation(mrg[i], mθg[j], mζg[k], ϕ[ind, :, :, :, :], rg, θg, ζg)
end

#%%
ft_int_phi = fft(int_phi, [2, 3]);

potential_plot(ft_int_phi, mgrids, label_max=0.3)

contour_plot(ϕ, grids, ind=ind)
contour_plot(int_phi, mgrids)
