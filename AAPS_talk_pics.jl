
#this file is used to generate some pictures for AAPPS.
#this should probably have been outside MID....

using MID
using MIDViz
using Plots

using LaTeXStrings


##############################
#cylindrical continuum.

Nr = 100;

geo = GeoParamsT(R0=4.0)

prob = init_problem(q=Axel_q, geo=geo)

rgrid = MID.ContGridDataT(N=Nr)
θgrid = asm_grid(N=2, start=2)
ζgrid = asm_grid(N=1, start=-2)

grids = init_grids(rgrid, θgrid, ζgrid)

evals_cont = cyl_cont(prob, grids);

continuum_plot(evals_cont)#, savefile="aapps_pics/cyl_cont.png")

############################################
#toroidal continuum. -> done a bit manually for clear labelling.

Nr = 100;

geo = GeoParamsT(R0=4.0)

prob = init_problem(q=Axel_q, geo=geo)

rgrid = MID.ContGridDataT(N=Nr)
θgrid = asm_grid(N=2, start=2)
ζgrid = asm_grid(N=1, start=-2)

grids = init_grids(rgrid, θgrid, ζgrid)




ω_cont = continuum(prob, grids);
#continuum_plot(evals, savefile="aapps_pics/")


rgrid_cont = MID.inst_grid(grids.r);

#guess -> good one lol.
split = 0.6
rgrid1 = rgrid[rgrid .< split]
rgrid2 = rgrid[rgrid .>= split]

ω_cont1 = ω_cont[rgrid .< split, :, :]
ω_cont2 = ω_cont[rgrid .>= split, :, :]

#niice
ωm2 = vcat(ω_cont1[:, 1, 1], ω_cont2[:, 2, 1])
ωm3 = vcat(ω_cont1[:, 2, 1], ω_cont2[:, 1, 1])

#groovy!
p = scatter(xlabel=L"r", ylabel=L"\frac{\omega  R_0}{v_A}", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10, ylimits=(-0.05, 1.05))

scatter!(rgrid, ωm2, label="(2, -2)")
scatter!(rgrid, ωm3, label="(3, -2)")
#scatter!(rgrid1, ω_cont1[:, 2, 1])

savefig(p, "aapps_pics/tor_cont.png")
#scatter!(rgrid2, ω_cont2[:, 1, 1])
#scatter!(rgrid2, ω_cont2[:, 2, 1])


#labels here are no good.
#continuum_plot(ωcont, grids)

#now we want to add the tae in!

#so we split the values up into 4 section.


Nr = 100;

geo = GeoParamsT(R0=4.0)

prob = init_problem(q=Axel_q, geo=geo)

rgrid = rfem_grid(N=Nr)
θgrid = asm_grid(N=2, start=2)
ζgrid = asm_grid(N=1, start=-2)

grids = init_grids(rgrid, θgrid, ζgrid)

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);

continuum_plot(evals)

tae_ind = find_ind(evals, 0.4)

potential_plot(ϕft, grids, 198, savefile="aapps_pics/tae_mode_structure.png")

rtae = rgrid_cont[argmax(abs.(ϕft[tae_ind, :, 1, 1]))]

ωtae = evals.ω[tae_ind]

#groovy!
p = scatter(xlabel=L"r", ylabel=L"\frac{\omega  R_0}{v_A}", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10, ylimits=(-0.05, 1.05))

scatter!(rgrid_cont, ωm2, label="(2, -2)")
scatter!(rgrid_cont, ωm3, label="(3, -2)")

scatter!([rtae], real.([ωtae]), label="TAE")

savefig("aapps_pics/tor_cont_tae.png")


##########################
#continuum mode.


Nr = 100;

geo = GeoParamsT(R0=4.0)

prob = init_problem(q=Axel_q, geo=geo, met=cylindrical_metric!)

rgrid = rfem_grid(N=Nr)
θgrid = asm_grid(N=2, start=2)
ζgrid = asm_grid(N=1, start=-2)

grids = init_grids(rgrid, θgrid, ζgrid)

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);

continuum_plot(evals)

mode_ind = find_ind(evals, 0.8)

display(evals.ω[mode_ind])
display(evals.r[mode_ind])

potential_plot(-1 .* ϕft, grids, mode_ind, savefile="aapps_pics/continuum_mode.png")



#############################
#toroidal correction continuum mode.


Nr = 100;

geo = GeoParamsT(R0=4.0)

prob = init_problem(q=Axel_q, geo=geo)

rgrid = rfem_grid(N=Nr)
θgrid = asm_grid(N=2, start=2)
ζgrid = asm_grid(N=1, start=-2)

grids = init_grids(rgrid, θgrid, ζgrid)

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);

continuum_plot(evals)

mode_ind = find_ind(evals, 0.8)

potential_plot(ϕft, grids, mode_ind, savefile="aapps_pics/continuum_mode_tor_correction.png")



####################################
#B mag picture
#given how fkn difficult this is, we may want to implement some square to circle mapping. Especially if we have Δ ≠ 0.
#and we would also like to show the coordinates.

function mod_B(r, θ)
    B0 = 1
    R0 = 4
    return B0 * (1 - r/R0 * cos(θ))

end

rB = LinRange(0, 1, 100)
θB = LinRange(0, 2π, 100)

magB = zeros(100, 100)

for (i, r) in enumerate(rB), (j, θ) in enumerate(θB)

    magB[i, j] = mod_B(r, θ)
end


RB = @. 4 + rB * cos(θB)

ZB = @. rB * sin(θB)

contourf(θB, rB, magB)
scatter(ZB, RB, magB)


heatmap(θB, rB,  mod_B.(rB, θB'), proj=:polar, axis=false, dpi=600)#,savefig="aapps_pics/magB_torus.png")
#contourf(θB, rB,  mod_B.(rB, θB')', proj=:polar)
savefig("aapps_pics/magB_torus.png")



############################################

#island coordinates.
#see restricted mapping.jl








###############################################

#island continuum -> also created a combined plot with Gadi datae

#island from (2, -1) case on Gadi.
isl = IslandT(2, -1, 0.00015625000000000003, 2.0, 2.0, 0.5, 0.05)

geo = GeoParamsT(R0=10.0)

#start with this???
pmd = asm_grid(start=-12, N=26, incr=1)
tmd = asm_grid(start=-2, N=5, incr=1)

κlist = LinRange(0.000001, 0.999, 100)

χlist = @. -(2*isl.A*κlist - isl.A)

ω2list = island_continuum(χlist, pmd, tmd, geo, isl, 0);


scatter(κlist, sqrt.(abs.(ω2list .* geo.R0^2)), legend=false, markersize=0.1)#, ylimits=(0.3, 0.5))#.2, 0.6))


###################################################
#unperturbed continuum.

Nr = 1000;

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=island_mode_21, geo=geo, met=cylindrical_metric!)

rgrid = MID.ContGridDataT(N=Nr)
θgrid = asm_grid(start=-1, N=9)
ζgrid = asm_grid(start=-4, N=5)

grids = init_grids(rgrid, θgrid, ζgrid)

evals_cont = cyl_cont(prob, grids);



p = scatter(evals_cont.r, real.(evals_cont.ω), group=evals_cont.modelabs, xlabel=L"r", ylabel=L"\frac{\omega  R_0}{v_A}", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10, ylimits=(-0.005, 0.08), legend=false, xlimits=(0.45, 0.55))

display(p)

savefig("aapps_pics/unperturbed_cont.png")

#continuum_plot(evals_cont, ymax=0.08, ymin=-0.005)#, savefile="aapps_pics/cyl_cont.png")