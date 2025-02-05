
#Want q-profile that is pretty flat, and gives low freq's to best match island
#want 2, 1 island at r=0.5, so q(0.5)=2
#use density to create v small bump to give gae.
#so q should be flat ish

#ok, so this combo gives a (2, -1) gae at ω=0.0075, peaking at ~0.68.
#think this will make a pretty good setup.

using MID
using Plots; plotlyjs()

function nq(r)
    a = 1.98
    b = 0.08

    q = a + b*r^2
    dq = 2 * b * r

    return q, dq
end

#this is perhaps going to be more complicated unfort.
#so this is not the go, rather we want to use something like gae_isl_dens, but maybe want to modify κ, however the q-profile is probably a pretty good start.
function ndens(r)

    p = 1
    κ = -2

    return ((1-p*r^2)^κ)

    """
    left = 0.5
    right = 0.7
    scale = 10

    if r<left
        return 1.0
    elseif r>right
        return 1.0
    else
        return scale*r^2 -scale * (left+right)*r + 1 + left*right*scale
    end
    """
end

rplot = LinRange(0, 1, 100)

plot(rplot, ndens.(rplot))

#plot(rplot, gae_isl_dens.(rplot))

plot(rplot, axel_dens.(rplot))

plot(rplot, nq.(rplot))

#with individual parts.
#different ways of doing it are basically identical.
#10.992975 seconds (15.38 M allocations: 1.195 GiB, 4.40% gc time)

#start very small, matrix scales much more extremly
Nr = 500;
Nθ = 8

geo = GeoParamsT(R0=1000.0)

#isl = IslandT(A=0.0e-4, m0=3, n0=2);

prob = init_problem(q=nq, geo=geo, dens=ndens)#, isl=isl); 

rgrid = rfem_grid(N=Nr)
#θgrid = afem_grid(N=Nθ, pf=3)
θgrid = asm_grid(start=0, N=3)
ζgrid = asm_grid(start=-1, N=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)

#cont_grids = init_grids(Nr, θgrid, ζgrid)

#ω_cont = continuum(prob, cont_grids)

#plot_continuum(ω_cont, cont_grids)

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=false, target_freq=0.008, nev=300);

#so island is ω ~ 0.05 in Axel's case
continuum_plot(evals, ymax=0.1)
tae_ind = find_ind(evals, 0.0075775)

potential_plot(ϕft, grids, tae_ind)

rgrid_cont = MID.ContGridDataT(N=Nr)
cont_grids = init_grids(rgrid_cont, θgrid, ζgrid)

evals_cont = cyl_cont(prob, cont_grids)

continuum_plot(evals_cont, ymax=0.1)

#perhaps with more res?

