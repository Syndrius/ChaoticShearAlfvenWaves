#we may also need to adjust slightly the chaotic cases to make the region they explore a bit
#more consistent!

#we probabbly will need to change the shape of the two irrational surfaces!
#neat copies of the poincare plots for each case.
using MID
using MIDCantori
using MIDViz
using Plots
using LaTeXStrings
using JLD2
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
#the two irrationals are v different now!
conv_num1 = 9 #10 is a big jump, but not impossible, perhaps for the final run?
conv_num2 = 11 #14 causes us to run out of mem, takes fkn ages.#probably want this to be 15 for actual results, finding the all the orbits takes ~20 mins, smaller number for changing the plot vars.
#stuart uses 22nd convergents, for that we will need Gadi. -> his case is also just above the critical threshold!
ylims = (0.5, 2/3)
tfs = 16
dpi = 600
msize = 0.7
msize1 = 1.7
msize2 = 1.7
#1 is for (6, 4) branch
c1 = cur_colors[2]
c2 = cur_colors[3]
#c1 =:black
#c2 =:black

shape1 = :o
shape2 = :o
#%%
#extracted from Gadi, using the 16th convergent.
#for ind2, this was (6249, 4351)
#this is obviously a stupid way of doing this, but the data was split across two runs.
x0s = zeros(2, 8)
#was done with k0, k08, k10, k11, k12, k13, k15, k17
#may need k05 as well!
x0s[:, 1] = [0.6074571931509041, 0.0]
x0s[:, 2] = [0.602560800274648, 2.2009024459434352e-9]
x0s[:, 3] = [0.6014998588000774, 1.2402710810156599e-10]
x0s[:, 4] = [0.6009927005823444, -1.1570874341706402e-9]
x0s[:, 5] = [0.6004988361156988, -1.0902676386923584e-9]
x0s[:, 6] = [0.6000166548016622, -1.5287005060937134e-10]
x0s[:, 7] = [0.5990874285426273, 1.742571748958684e-11]
x0s[:, 8] = [0.5982316368030925, 3.021340507842228e-12]
#%%
MIDCantori.NumberTheory.continued_fraction(1.4362215818857254, 15)
#probably not perf but good enough!
ir1 = [1, 2, 1, 2, 2, 1, 5, 1, 5, 1, 5, 1, 5, 1, 5, 1, 5, 1, 5, 1, 5, 1, 5, 1] #not very noble!
ir2 = [1, 2, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
#conv1 = MIDCantori.NumberTheory.convergents(ir1)
#a1 = conv1[conv_num][2]
#b1 = conv1[conv_num][1]
conv1 = MIDCantori.NumberTheory.convergents(ir1)
conv2 = MIDCantori.NumberTheory.convergents(ir2)
#think the order of this depends on how we define our irrational number
#keeping up the incosistency with a/b!
a1 = conv1[conv_num1][1]
b1 = conv1[conv_num1][2]
a2 = conv2[conv_num2][1]
b2 = conv2[conv_num2][2]

Ntraj = 151
ψvals = collect(LinRange(0.35, 0.8, Ntraj));
θvals = zeros(Ntraj);
k0 = 0.0
k05 = 0.0005 #shit example
k08 = 0.0008
k10 = 0.0010
k11 = 0.0011
k12 = 0.0012
k13 = 0.0013
k15 = 0.0015
k17 = 0.0017
#%%
k = k0
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)

r01, θ01 = periodic_orbit(a1, b1, 4, k)
r02, θ02 = periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 5000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#x1, z1 = poincare_plot(prob, 2000, [r02[1]], [θ02[1]]);
#x2, z2 = poincare_plot(prob, 2000, [x0s[1, 1]], [x0s[2, 1]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k0_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
#scatter!(x2[1, :], z2[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k0_poincare.png")
#%%
#k05 is a shit example.
k = k08
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)

#r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r01, θ01 = periodic_orbit(a1, b1, 4, k)
r02, θ02 = periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 5000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#x1, z1 = poincare_plot(prob, 2000, [r02[1]], [θ02[1]]);
#x2, z2 = poincare_plot(prob, 2000, [x0s[1, 2]], [x0s[2, 2]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k08_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
#scatter!(x2[1, :], z2[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k08_poincare.png")

#%%
#this case is probably a bit boring!
k = k10
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)

#r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r01, θ01 = periodic_orbit(a1, b1, 4, k)
r02, θ02 = periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 5000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#x1, z1 = poincare_plot(prob, 2000, [r02[1]], [θ02[1]]);
#x2, z2 = poincare_plot(prob, 2000, [x0s[1, 3]], [x0s[2, 3]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k10_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
#scatter!(x2[1, :], z2[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k10_poincare.png")
#%%
#could be a tricky one to explain, as the first flux surface is broken (visble!)
#but it does not yet go all over the place!
k = k11
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)

#r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r01, θ01 = periodic_orbit(a1, b1, 4, k)
r02, θ02 = periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
#x2, z2 = poincare_plot(prob, 2000, [x0s[1, 4]], [x0s[2, 4]]);
x1, z1 = poincare_plot(prob, 20000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#x1, z1 = poincare_plot(prob, 2000, [r02[1]], [θ02[1]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k11_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
#scatter!(x2[1, :], z2[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k11_poincare.png")
#%%
k = k12
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)

#r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r01, θ01 = periodic_orbit(a1, b1, 4, k)
r02, θ02 = periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 20000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#x1, z1 = poincare_plot(prob, 2000, [r02[1]], [θ02[1]]);
#x2, z2 = poincare_plot(prob, 2000, [x0s[1, 5]], [x0s[2, 5]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k12_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
#scatter!(x2[1, :], z2[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k12_poincare.png")

#%%
k = k13
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)

#r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r01, θ01 = periodic_orbit(a1, b1, 4, k) #now way above critical value, starts to become difficult. #can probably just pick a roughly close value!
r02, θ02 = periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 20000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#x1, z1 = poincare_plot(prob, 2000, [r02[1]], [θ02[1]]);
#x2, z2 = poincare_plot(prob, 2000, [x0s[1, 6]], [x0s[2, 6]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k13_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
#scatter!(x2[1, :], z2[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k13_poincare.png")

#%%
k = k15
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)
#r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r01, θ01 = periodic_orbit(a1, b1, 4, k)
r02, θ02 = periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 20000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#not super sure what is best for this one.
#x1, z1 = poincare_plot(prob, 20000, [r02[1]], [θ02[1]]);
#x2, z2 = poincare_plot(prob, 20000, [x0s[1, 7]], [x0s[2, 7]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k15_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
#scatter!(x2[1, :], z2[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k15_poincare.png")
#%%
k = k17
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)
#r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
r01, θ01 = periodic_orbit(a1, b1, 4, k)
r02, θ02 = periodic_orbit(a2, b2, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);
x1, z1 = poincare_plot(prob, 20000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#ok so we just have to run this for longer!
#now we see that the orbit can go basically anywhere.
#x1, z1 = poincare_plot(prob, 10000, [r02[1]], [θ02[1]]);
#x2, z2 = poincare_plot(prob, 10000, [x0s[1, 8]], [x0s[2, 8]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/k17_no_ir_poincare.png")
scatter!(x1[1, :], z1[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
#scatter!(x2[1, :], z2[1, :], markersize=msize1, color=c1, label=L"\iota_1", markershape=shape1)
scatter!(x1[2, :], z1[2, :], markersize=msize2, color=c2, label=L"\iota_2", markershape=shape2)
savefig("~/phd/MID/QFMPaper/results/k17_poincare.png")
#%%
