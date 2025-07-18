
using MID
using MIDViz
using JLD2
using Plots; plotlyjs()
using Plots; gr()
using Statistics
#%%
#based on umbrella, last KAM surface breaks at ~0.0014
#bit worried about the huge islands in the middle.
k_min = 0.0013 #This case should have a few solid surfaces.

k_mid = 0.00145 #no solid KAM surfaces, but clearly still structure.

k_max = 0.0017 #above this the qfm surfaces inside the chaotic region get cooked, this may still be too much though!

k = k_min
k1 = 0.00115
k15 = 0.00125
k2 = 0.0013
k3 = 0.00145
k4 = 0.0017
k0 = 0.0005
k = k2

geo = init_geo(R0=1.0)
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)

isls = MID.Geometry.IslandT[isl1, isl2]

prob = init_problem(geo=geo, q=cantori_q, isls=isls, met=:cylinder, type=:flux)

M = 32
N = 8

#%%

Ntraj = 300;
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
ψlist = collect(LinRange(0.4, 0.8, Ntraj));
Nlaps = 500;

#would be good if we could plot the poincare with like a gradient of colours to see where they end up.
poincare_plot(prob, Nlaps, ψlist, zeros(Ntraj), ylimits=(0.5, 0.67))

#%%

total_rats = [lowest_rationals(11, prob.q(0.0)[1], prob.q(1.0)[1]);
              lowest_rationals(21, prob.q(0.96)[1], prob.q(1.0)[1]);
              lowest_rationals(25, prob.q(0.0)[1], prob.q(0.15)[1]);
              [(52, 51), (44, 43), (31, 30)]];

total_gl = surface_guess(total_rats, prob.q)
total_surfs = construct_surfaces(total_rats, total_gl, prob, M=M, N=N);
plot_surfs(total_surfs)
#%%
#extras inside the chaotic region
#note that some of these are not there!
extra_rats1 = [(17, 12), (18, 13), (20, 13)]
extra_rats2 = [(19, 13), (19, 14), (17, 13)]
egl1 = surface_guess(extra_rats1, prob.q)
egl2 = surface_guess(extra_rats2, prob.q)
esurfs1 = construct_surfaces(extra_rats1, egl1, prob, M=M, N=N);
esurfs2 = construct_surfaces(extra_rats2, egl2, prob, M=M, N=N);
extra_rats3 = []
plot_surfs(esurfs1)

#%%
#divying up the surfaces into increments so we can show how B and J change with more surfaces.
rats1 = [lowest_rationals(4, prob.q(0.0)[1], prob.q(1.0)[1]); (52, 51)]
inds1 = [i for i in 1:length(total_rats) if total_rats[i] in rats1]
surfs1 = total_surfs[inds1];
rats2 = [lowest_rationals(6, prob.q(0.0)[1], prob.q(1.0)[1]); (52, 51)]
inds2 = [i for i in 1:length(total_rats) if total_rats[i] in rats2];
surfs2 = total_surfs[inds2];
rats3 = [lowest_rationals(8, prob.q(0.0)[1], prob.q(1.0)[1]); (52, 51)]
inds3 = [i for i in 1:length(total_rats) if total_rats[i] in rats3];
surfs3 = total_surfs[inds3];
rats4 = [lowest_rationals(10, prob.q(0.0)[1], prob.q(1.0)[1]); (52, 51)]
inds4 = [i for i in 1:length(total_rats) if total_rats[i] in rats4];
surfs4 = total_surfs[inds4];
#rats5 = total_rats;
rats5 = [(12, 11), (13, 11), (14, 11), (18, 11), (19, 11), (20, 11), (21, 11)]
inds5 = [i for i in 1:length(total_rats) if total_rats[i] in rats5];
surfs5 = total_surfs[inds5];
rats6 = [(15, 11), (16, 11), (17, 11)]
inds6 = [i for i in 1:length(total_rats) if total_rats[i] in rats6];
surfs6 = total_surfs[inds6];

#rats6 = lowest_rationals(21, prob.q(0.96)[1], prob.q(1.0)[1])
#rats7 = lowest_rationals(25, prob.q(0.0)[1], prob.q(0.15)[1])
#lowest_rationals(11, prob.q(0.0)[1], prob.q(1.0)[1])

#%%
using LaTeXStrings
using Plots; gr()
#ok this actually matches our finding really well
#surfl1 = surfs1;
#surfl2 = vcat(surfs1, surfs2);
#surfl3 = vcat(surfl2, surfs3);
#surfl4 = vcat(surfl3, surfs4);
#surfl5 = vcat(surfl4, surfs5);

#this time each iteration already includes the last
surfl1 = surfs1;
surfl2 = surfs2;
surfl3 = surfs3;
surfl4 = surfs4;
surfl5 = vcat(surfs4, surfs5);
surfl6 = vcat(surfl5, surfs6);
surfl7 = vcat(esurfs1, surfl6);
surfl8 = vcat(esurfs2, surfl7);
#surfl6 = vcat(surfl5, surfs7);
#surfl7 = vcat(surfl6, surfs6);
#surfl8 = vcat(surfl7, surfs8);
surfs_list = [surfl1, surfl2, surfl3, surfl4, surfl5, surfl6, surfl7, surfl8];

#we may want to change the first set of surfaces, as they artificially look better I think.
#note that Bs is actually ~10 times smaller at 6 vs 1.
#note that Helander et al uses ~60 surfaces, so perhaps we will want more?
nsurfs = [length(surfs) for surfs in surfs_list]
nsurfs
#%%
rgrid_jac = init_grid(type=:rf, N = 100, start=0.15, stop=0.95)
θgrid_jac = init_grid(type=:af, N = 20) 
ζgrid_jac = init_grid(type=:af, N = 4)
grids_jac = init_grids(rgrid_jac, θgrid_jac, ζgrid_jac)
#B, jac, djac = compute_jac(prob, grids_jac, surfs3);

#%%

B2mean = zeros(length(surfs_list));
djacdrmean = zeros(length(surfs_list));
djacdθmean = zeros(length(surfs_list));

for i in 1:length(surfs_list)
    display(i)
    curr_surfs = surfs_list[i]

    B, jac, djac = compute_jac(prob, grids_jac, curr_surfs);
    display("where prob")

    B2mean[i] = mean(B[1, :, :, :] .^ 2)
    djacdrmean[i] = mean(abs.(djac[1, :, :, :]))
    djacdθmean[i] = mean(abs.(djac[2, :, :, :]))
end
#%%
plot_font = "Computer Modern"
default(fontfamily=plot_font, grid=false, framestyle=:semi, palette=:tol_bright)
cur_colors = palette(:tol_bright);
lfs = 18
xfs = 20
tfs = 16
ylims = (0.0, 1)
xlims = (15, 45)
plot(nsurfs, B2mean, palette=:tol_bright, axis=:left, xlimits=xlims, label=L"(B^s)^2", color=cur_colors[1], legendfontsize=lfs, xlabel="Number of QFM Surfaces", xguidefontsize=xfs, xtickfontsize=tfs, ytickfontsize=tfs-3, dpi=1200, legend=:topleft, margin=5Plots.mm, y_foreground_color_border=cur_colors[1], y_foreground_color_axis=cur_colors[1], y_foreground_color_text=cur_colors[1], yaxis=:log)

#plot!(nsurfs, djacdrmean, axis=:right, palette=:tol_bright, label=L"\partial\mathcal{J}/\partial s")
plot!(twinx(), nsurfs, djacdrmean, palette=:tol_bright, xlimits=xlims, axis=:right, linestyle=:dot,linewidth=3, legend=:topright, label=L"|\partial\mathcal{J}/\partial s|", color=cur_colors[2], legendfontsize=lfs, xguidefontsize=xfs, tickfontsize=tfs, dpi=1200, y_foreground_color_border=cur_colors[2], y_foreground_color_axis=cur_colors[2], y_foreground_color_text=cur_colors[2], markerstrokewidth=10)#, yticklabelcolor=cur_colors[2])

plot!(twinx(), nsurfs .+ 100, djacdrmean, xlimits=xlims,  palette=:tol_bright, axis=:right, linestyle=:dot,linewidth=1.5, legend=:topright, label=L"|\partial\mathcal{J}/\partial s|", color=cur_colors[2], legendfontsize=lfs, xguidefontsize=xfs, tickfontsize=tfs, dpi=1200, y_foreground_color_border=cur_colors[2], y_foreground_color_axis=cur_colors[2], y_foreground_color_text=cur_colors[2], markerstrokewidth=10)#, yticklabelcolor=cur_colors[2])

savefig("/Users/matt/phd/MID/QFMPaper/results/number_of_surfaces_log.png")
#%%
#plot(nsurfs, djacdθmean)
#why the fek didn't we save the surf_list!!
#may want to redo this as a log plot.
save_object("/Users/matt/phd/QFMPaper/nsurfs.jld2", nsurfs)
save_object("/Users/matt/phd/QFMPaper/surfs_list.jld2", surfs_list)
save_object("/Users/matt/phd/QFMPaper/B2mean.jld2", B2mean)
save_object("/Users/matt/phd/QFMPaper/djacdrmean.jld2", djacdrmean)
good_rationals = vcat(rats1, rats2, rats3, rats4, rats5, rats6, rats7, rats8);
save_object("/Users/matt/phd/QFMPaper/rationals.jld2", good_rationasl)
#%%
nsurfs = load_object("/Users/matt/phd/QFMPaper/nsurfs.jld2");
surfs_list = load_object("/Users/matt/phd/QFMPaper/surfs_list.jld2");
B2mean = load_object("/Users/matt/phd/QFMPaper/B2mean.jld2");
djacdrmean = load_object("/Users/matt/phd/QFMPaper/djacdrmean.jld2");
#good_rationals = vcat(rats1, rats2, rats3, rats4, rats5, rats6, rats7, rats8);
good_rationals = load_object("/Users/matt/phd/QFMPaper/rationals.jld2");
#%%
nsurfs = load_object("/Users/matt/phd/QFMPaper/nsurfs.jld2")
surfs_list = load_object("/Users/matt/phd/QFMPaper/surfs_list.jld2");
B2mean = load_object("/Users/matt/phd/QFMPaper/B2mean.jld2")
djacdrmean = load_object("/Users/matt/phd/QFMPaper/djacdrmean.jld2")
#good_rationals = vcat(rats1, rats2, rats3, rats4, rats5, rats6, rats7, rats8);
#save_object("/Users/matt/phd/QFMPaper/rationals.jld2", good_rationasl)
#%%
