#neat copies of the qfm plots, including generic poincare plot
#overlay of surfs
#and qfm poincare plot.
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
conv_num = 13 #14 causes us to run out of mem, takes fkn ages.#probably want this to be 15 for actual results, finding the all the orbits takes ~20 mins, smaller number for changing the plot vars.
#stuart uses 22nd convergents, for that we will need Gadi. -> his case is also just above the critical threshold!
ylims = (0.45, 0.75)
tfs = 16
dpi = 300
msize = 0.7
msize1 = 1.7
msize2 = 1.7
c1 = cur_colors[1]
c2 = cur_colors[2]
#c1 =:black
#c2 =:black
ψvals = unique([collect(LinRange(0.3, 0.55, 30)) ; collect(LinRange(0.55, 0.62, 50)) ; collect(LinRange(0.62, 0.8, 20))])
θvals = zeros(length(ψvals))

shape1 = :o
shape2 = :o
#%%
#surface overlay plot
#do k13 case, as it has chaos, but also normal flux surfaces, probably best case for qfm poincare.
#note k2 is k13!
surfs = load_object("/Users/matt/phd/MID/data/surfaces/qfm/k2_surfs.jld2");
k = 0.0013
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geometry(:cyl, R0=1.0)
fields = init_fields(:ψ, q=cantori_q, isls=isls)
prob = init_problem(geometry=geo, fields=fields)
#r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
x, z = poincare_plot(prob, 1000, ψvals, θvals);

#x1, z1 = poincare_plot(prob, 2000, [r01[1], r02[1]], [θ01[1], θ02[1]]);
#%%
scatter(x, z, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/base_poincare.png")
plot_surfs(surfs, overlay=true, color=cur_colors[2], linewidth=3.0, ylimits=(0.45, 0.72))
savefig("~/phd/MID/QFMPaper/results/overlay_poincare.png")
#%%
xq, zq = poincare_plot(prob, 1000, ψvals, θvals, surfs=surfs);
#%%
scatter(xq, zq, ylimits=(0.45, 0.72), markersize=msize, legend=false, color=:black, xlabel=L"\vartheta", ylabel=L"s", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/qfm_poincare.png")
#%%
#unpert poincare
ψ0vals = LinRange(0.3, 0.8, 100)
θ0vals = zeros(length(ψ0vals))
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, type=:flux)
#r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
x0, z0 = poincare_plot(prob, 1000, ψ0vals, θ0vals);
#%%
scatter(x0, z0, ylimits=ylims, markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/unpert_poincare.png")
#%%
scatter(x, z, ylimits=(0.55, 0.63), markersize=msize, legend=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
#scatter(x, z, ylimits=ylims, markersize=msize, label=false, color=:black, xlabel=L"\theta", ylabel=L"\psi", xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, dpi=dpi)
savefig("~/phd/MID/QFMPaper/results/base_poincare_zoom.png")
