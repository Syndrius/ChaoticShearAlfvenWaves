
#neat copies of the poincare plots for each case.
using MID
using MIDCantori
using MIDViz
using Plots
using LaTeXStrings
#%%
plot_font = "Computer Modern"
#default(fontfamily=plot_font, grid=false, framestyle=:semi, palette=:tol_bright)
default(fontfamily=plot_font, grid=false, palette=:tol_bright)
cur_colors = palette(:tol_bright);
#lfs = 16
xfs = 28
yfs = 26
ytfs = 14
xtfs = 16
ylab = L"\psi"
xlab = L"\theta"
conv_num = 15 #probably want this to be 15 for actual results, finding the all the orbits takes ~20 mins, smaller number for changing the plot vars.
ylims = (0.5, 2/3)
tfs = 16
dpi = 1200
msize = 0.7
msize1 = 1.7
msize2 = 1.7
c1 = cur_colors[1]
c2 = cur_colors[2]
#c1 =:black
#c2 =:black

shape1 = :o
shape2 = :o
#%%
ir1 = [0, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
ir2 = [0, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
conv1 = MIDCantori.NumberTheory.convergents(ir1)
a1 = conv1[conv_num][2]
b1 = conv1[conv_num][1]
conv2 = MIDCantori.NumberTheory.convergents(ir2)
a2 = conv2[conv_num][2]
b2 = conv2[conv_num][1]

Ntraj = 151
ψvals = collect(LinRange(0.35, 0.8, Ntraj));
θvals = zeros(Ntraj);
#%%
k = 0.0005
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)

r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r02, θ02 = MIDCantori.QFM.anal_periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 2000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k0_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k0_poincare.png")

#%%
k = 0.00115
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)

r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r02, θ02 = MIDCantori.QFM.anal_periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 2000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k1_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k1_poincare.png")

#%%
k = 0.00125
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)

r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r02, θ02 = MIDCantori.QFM.anal_periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 2000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k15_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k15_poincare.png")

#%%
k = 0.0013
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)
r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r02, θ02 = MIDCantori.QFM.anal_periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 2000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k2_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k2_poincare.png")
#%%
k = 0.00145
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)
r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r02, θ02 = MIDCantori.QFM.anal_periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 2000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k3_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k3_poincare.png")
#%%
k = 0.0017
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)
r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r02, θ02 = MIDCantori.QFM.anal_periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 2000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k4_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k4_poincare.png")
#%%
k = 0.0
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)
r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r02, θ02 = MIDCantori.QFM.anal_periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 2000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k0_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k0_poincare.png")
