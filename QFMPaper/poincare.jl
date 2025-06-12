#here we generate the nice poincae plots for puplication.

using MID
using MIDViz
#using LaTeXStrings
using JLD2
using Plots; gr()
#%%
function cyl_qfm_q(ψ::Float64)

    return 1.0 + ψ, 1.0
end
#%%
geo = init_geo(R0=10.0)
isl1 = init_island(m0=7, n0=-5, A=0.006/7)
isl2 = init_island(m0=3, n0=-2, A=0.004/3)
isl3 = init_island(m0=8, n0=-5, A=0.006/8)

isls = [isl1, isl2, isl3]

prob = init_problem(geo=geo, q=cyl_qfm_q, isls=isls, met=:cylinder)
#%%
Ntraj = 100
rlist = collect(LinRange(0.2, 0.8, Ntraj));
#rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 800;

x, z = poincare_plot(prob, Nlaps, Ntraj, rlist)#, color=:black)

#%%

surfs = load_object("/Users/matt/phd/low_shear_surfs.jld2");
length(surfs)
#plot_surfs(surfs)
#%%

Ntraj = 80
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
#changing this range, there is still the weird 3 island at ~0.68
#that is exactly where the (3,, 2) island would be, guess this is an example of the islands combining?
rlist = collect(LinRange(0.5, 0.95, Ntraj));
Nlaps = 500;

#this looks awful? Are we doing something wrong here?
#odd that the benchmark case seems to be working perfectly
#does kind of look like Stuart never shows the qfm surfaces being able to recreate a straight poincare section
#at least ot in the chaotic region.
#his does look much much nicer though.
#seems especially worrying that the new poincare section is cooked away from the chaotic region.
#perhaps this signifies deeper problemos?
#the region we were going to show is pretty much fine, however, the extra bit has a weird extra island.
xs, zs = poincare_plot(prob, Nlaps, Ntraj, rlist, surfs=surfs)

#%%
surface_guess([(3, 2)], prob.q)


plot_font = "Computer Modern"
default(fontfamily=plot_font, grid=false, framestyle=:semi, palette=:tol_bright)
cur_colors = palette(:tol_bright);
lfs = 18
xfs = 22
yfs = 22
tfs = 16
ylims = (0.7, 0.95)
MIDViz.plot_poincare(x, z, color=:black, xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=tfs, ytickfontsize=tfs-2, dpi=1200, ylimits=ylims)
#not really sure if this actually wants the ticklabels tbh.
savefig("/Users/matt/phd/QFMPaper/poincare.png")
#%%
MIDViz.plot_poincare(xs, zs, xlabel=L"ϑ", ylabel=L"s", color=:black, xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=tfs, ytickfontsize=tfs-2, dpi=1200, ylimits=ylims)

savefig("/Users/matt/phd/QFMPaper/poincare_qfm.png")
#%%
MIDViz.plot_poincare(x, z, color=:black, xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=tfs, ytickfontsize=tfs-2, dpi=1200, ylimits=ylims)
MIDViz.Plotting.overlay_surfs(surfs, color=cur_colors[2], linewidth=2.5)
savefig("/Users/matt/phd/QFMPaper/poincare_overlay.png")
#%%
plot_surfs(surfs, ylimits=(0.7, 0.95), legend=false, color=:black, dpi=1200, xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=tfs, ytickfontsize=tfs-2)
savefig("/Users/matt/phd/QFMPaper/surfs.png")

