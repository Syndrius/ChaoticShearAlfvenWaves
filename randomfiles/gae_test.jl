
#lets see if we can have a look at a gae!
#found one, should probably incorporate this properly once we understand this!

using MID
using Plots; plotlyjs()

#about as big as direct solving can handle
#for this case, fss is much better, 
Nr = 100
Nθ = 10


geo = GeoParamsT(R0 = 1000)
rgrid = rfem_grid(N=60, start=0.0, stop=0.999)
θgrid = asm_grid(start=-2, N=5)
#θgrid = afem_grid(N=10, pf=1)
ζgrid = asm_grid(start=0, N=2)

grids = init_grids(rgrid, θgrid, ζgrid)

function isl_q(r::Float64)
    q0 = 2/1
    qp = 2 #chosen pretty arbitrarily based on vibes of continuum.
    qp = -8
    r0 = 0.5
    #ψ0 = 0.125
    #q = 1 / (1 / q0 - qp / (q0)^2 * (r^2/2-ψ0))
    #dq = 4*(q0)^2* qp * r / (2*(q0)-qp*r^2 + 2*qp*ψ0)^2
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (r^2-r0^2))
    dq = 4*qp*q0^2*r*r0 / (2*q0*r0 - qp * (r^2-r0^2))^2
    return q, dq
end

function gae_q(r)

    #paper says constant, that may be after the Geometry contribution though!

    q = 1/(r^2)
    dq = -2.0/r^3
    return q, dq
end

function test_gae_dens(r)
    κ = -1 #paper just says less than 0?
    p = 1
    q0 = 2/1
    qp = 2 #chosen pretty arbitrarily based on vibes of continuum.
    r0 = 0.5

    g0 = 0.8
    #p = qp / (2 * g0^2*qp - 2*q0*r0 - qp*r0^2)

    return ((1-p*r^2)^κ)
end

#prob = init_problem(q = gae_q, dens=gae_dens, geo=geo)
#isl = IslandT(m0=2, n0=-1, w=0.1)
prob = init_problem(q = gae_isl_q, geo=geo, met=cylindrical_metric!, dens=gae_isl_dens)#, isl=isl)
prob = init_problem(q = gae_q, geo=geo, met=cylindrical_metric!, dens=test_gae_dens)


evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true, target_freq = 0.35, nev=100);

continuum_plot(evals)#, ymax=10)


#so we are getting like 10 gae's... wild.
#this might work, problem being the 10 gae's we are finding is a bit annoying!

#pretty similar results with both methods, fss is significanlty faster though!
gae_ind = find_ind(evals, 0.3543694)
tae_ind = find_ind(ω, 0.49)

gae_ind = 8

1.5/ 0.0003


potential_plot(4000 .* ϕft, grids, gae_ind)

using LaTeXStrings
using Printf
using Plots; gr()

p = plot(xlabel=L"r", ylabel=L"\phi", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)


rgrid, θgrid, _ = inst_grids(grids)



for i in 1:grids.θ.N

    for j in 1:grids.ζ.N
    

        

        plot!(rgrid,  real.(ϕft[gae_ind, :, i, j]), label=false)

    end

end

plot!(rgrid, 4000 .* real.(ϕft[gae_ind, :, 2, 1]), label=@sprintf("(1, 0)"), color=1)

#display(p)

savefig("aapps_pics/unpertuebd_gae.png")



#what about the analytical continuum!

rplot = LinRange(0.0001, 1, 100)

qplot = getfield.(gae_q.(rplot), 1)

qplot_isl = getfield.(isl_q.(rplot), 1)

nplot = gae_dens.(rplot)




plot(rplot, @. 1/qplot / (nplot))
plot(rplot, @. 1/qplot_isl / sqrt(nplot))

plot(rplot, nplot)

plot(rplot, qplot, ylimits=(0, 10))
