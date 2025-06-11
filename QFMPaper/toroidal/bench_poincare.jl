#here we generate the nice poincae plots for puplication.

using MID
using MIDViz
using LaTeXStrings
using JLD2
using Plots; gr()
#%%

k1 = 0.02
geo = init_geo(R0=4.0)
isl1 = init_island(m0=3, n0=2, A=k1/3)

isls = [isl1]#, isl2, isl3]

prob = init_problem(geo=geo, q=low_shear_qfm_q, isls=isls)
#%%
Ntraj = 80
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

x, z = poincare_plot(prob, Nlaps, Ntraj, rlist)

#%%

surfs = load_object("/Users/matt/phd/bench_low_shear_surfs.jld2");
length(surfs)
#%%

Ntraj = 80
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

xs, zs = poincare_plot(prob, Nlaps, Ntraj, rlist, surfs=surfs)

#%%


plot_font = "Computer Modern"
default(fontfamily=plot_font, grid=false, framestyle=:semi, palette=:tol_bright)
cur_colors = palette(:tol_bright);
lfs = 18
xfs = 22
yfs = 22
tfs = 16
ylims = (0.3, 0.7)
MIDViz.plot_poincare(x, z, color=:black, xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=tfs, ytickfontsize=tfs-2, dpi=1200, ylimits=ylims)
#not really sure if this actually wants the ticklabels tbh.
#savefig("/Users/matt/phd/QFMPaper/bench_poincare.png")
#%%
MIDViz.plot_poincare(xs, zs, xlabel=L"ϑ", ylabel=L"s", color=:black, xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=tfs, ytickfontsize=tfs-2, dpi=1200, ylimits=ylims)

#savefig("/Users/matt/phd/QFMPaper/bench_poincare_qfm.png")
#%%
MIDViz.plot_poincare(x, z, color=:black, xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=tfs, ytickfontsize=tfs-2, dpi=1200, ylimits=ylims)
MIDViz.Plotting.overlay_surfs(surfs, color=cur_colors[2], linewidth=2.5)
savefig("/Users/matt/phd/QFMPaper/bench_poincare_overlay.png")
#%%
plot_surfs(surfs, ylimits=(0.3, 0.7), legend=false, color=:black, dpi=1200, xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=tfs, ytickfontsize=tfs-2)
#savefig("/Users/matt/phd/QFMPaper/bench_surfs.png")

