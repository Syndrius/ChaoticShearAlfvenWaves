#should change the name of this fkn file, actually contains nice plotting of the B and jac vs number of surfaces.

#want to find a q-profile where there is a tae between the islands,
#as this should give it the best chance of resolving nicely.

using MID
using MIDViz
#using Plots; plotlyjs()
using Plots; gr()
using Statistics
using LaTeXStrings
#%%

function q_prof(r::Float64)
    #chosen so (4, 3)=0.4, (3, 2)=0.6
    #sillier values are same but with flux surfaces shifted from 0, 1 to (0.05, 0.95)
    a = 6/5
    a = 539/450
    b = 5/6
    b = 8/9
    return a + b*r^2, 2 * b * r
end

#%%

rgrid = init_grid(type=:rf, N = 150, sep1=0.4, sep2=0.6, frac=0.5, start=0.05, stop=0.95)
θgrid = init_grid(type=:as, start=-1, N = 6)
ζgrid = init_grid(type=:as, start=-5, N = 5)

grids = init_grids(rgrid, θgrid, ζgrid)

#%%

geo = init_geo(R0=4.0)
k = 0.001
isl = init_island(m0=4, n0=-3, A = k/4)
isl2 = init_island(m0=3, n0=-2, A = k/2)
prob = init_problem(geo=geo, q=q_prof, isl=isl, isl2=isl2)
unprob = init_problem(geo=geo, q=q_prof)

#%%

solver = init_solver(nev=200, targets=[0.23, 0.3, 0.37], prob=unprob)

#%%

evals_norm, ϕ_norm, ϕft_norm = compute_spectrum(prob=unprob, grids=grids, solver=solver);
#%%

continuum_plot(evals_norm)
ind_norm = find_ind(evals_norm, 0.2944)

potential_plot(ϕft_norm, grids, ind_norm)

#%%

Ntraj = 100;
#rlist = collect(LinRange(0.45, 0.65, Ntraj));
rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

#poincare_plot(qfm_benchmark_q, slab_to_plot, Nlaps, Ntraj, 0, 0.0, 0.0, 0.0, R0, isl, MID.Structures.no_isl, rlist)
poincare_plot(prob, Nlaps, Ntraj, rlist,  prob.geo.R0)#, filename=save_dir * "original_poincare.png");
q_prof(1.0)
#%%
rationals = [(6, 5), (2, 1), (3, 2), (4, 3), (5, 3), (5, 4), (7, 4), (7, 5), (8, 5), (9, 5)]

qvals = [sqrt(i[1] / i[2] - 539/450) / sqrt(8/9) for i in rationals]
qlist = [i[1] for i in rationals]
plist = [i[2] for i in rationals]

surfs = construct_surfaces(plist, qlist, qvals, prob);

plot_surfs(surfs)

#%%
rationals2 = [(11, 9), (11, 6)]
qvals2 = [sqrt(i[1] / i[2] - 539/450) / sqrt(8/9) for i in rationals2]
qlist2 = [i[1] for i in rationals2]
plist2 = [i[2] for i in rationals2]
surfs2 = construct_surfaces(plist2, qlist2, qvals2, prob);
#significant improvment, even to the region far from these two surfaces.
plot_surfs(surfs2)
plot_surfs(vcat(surfs, surfs2))
#%%
rationals3 = [(9, 7), (11, 7)]
qvals3 = [sqrt(i[1] / i[2] - 539/450) / sqrt(8/9) for i in rationals3]
qlist3 = [i[1] for i in rationals3]
plist3 = [i[2] for i in rationals3]
surfs3 = construct_surfaces(plist3, qlist3, qvals3, prob);
#made minimal difference to spectrum, probbaly shows we need more in the chaotic region
plot_surfs(surfs3)
plot_surfs(vcat(surfs, surfs2, surfs3))
#%%
#add two in the chaotic region and one close to zero.
rationals4 = [(17, 14), (11, 8), (10, 7)]
qvals4 = [sqrt(i[1] / i[2] - 539/450) / sqrt(8/9) for i in rationals4]
qlist4 = [i[1] for i in rationals4]
plist4 = [i[2] for i in rationals4]
surfs4 = construct_surfaces(plist4, qlist4, qvals4, prob);
#TAE is actually somewhat visible now, and looks like a tae, shows the 4,3 island is a bit bananas
#also area surrounding tae is crazy
plot_surfs(surfs4)
plot_surfs(vcat(surfs, surfs2, surfs3, surfs4))
#%%
#highest jac, and problems seem to be occuring just above r=0.4 and just below r=0.6
#add surfaces there.
rationals5 = [(15, 11), (13, 9)]
qvals5 = [sqrt(i[1] / i[2] - 539/450) / sqrt(8/9) for i in rationals5]
qlist5 = [i[1] for i in rationals5]
plist5 = [i[2] for i in rationals5]
surfs5 = construct_surfaces(plist5, qlist5, qvals5, prob);
#improvment around tae, but looks ot be worse at r=0.55 ish. and perhaps at r=0.25 ish
plot_surfs(surfs5)
plot_surfs(vcat(surfs, surfs2, surfs3, surfs4, surfs5))
#%%
#highest jac, and problems seem to be occuring at r=0.4 and at r=0.55
#add surfaces there.
rationals6 = [(13, 10), (14, 9)]
qvals6 = [sqrt(i[1] / i[2] - 539/450) / sqrt(8/9) for i in rationals6]
qlist6 = [i[1] for i in rationals6]
plist6 = [i[2] for i in rationals6]
surfs6 = construct_surfaces(plist6, qlist6, qvals6, prob);
#improvment around tae, but looks ot be worse at r=0.55 ish. and perhaps at r=0.25 ish
plot_surfs(surfs6)
plot_surfs(vcat(surfs, surfs2, surfs3, surfs4, surfs5, surfs6))
#%%
#adding some more for low r, as there seems to be other issues there.
rationals7 = [(16, 13), (23, 19)]
qvals7 = [sqrt(i[1] / i[2] - 539/450) / sqrt(8/9) for i in rationals7]
qlist7 = [i[1] for i in rationals7]
plist7 = [i[2] for i in rationals7]
#these are very good surfaces
surfs7 = construct_surfaces(plist7, qlist7, qvals7, prob);
#improvment around tae, but looks ot be worse at r=0.55 ish. and perhaps at r=0.25 ish
plot_surfs(surfs7)
plot_surfs(vcat(surfs, surfs2, surfs3, surfs4, surfs5, surfs6, surfs7))
#%%
rationals8 = [(17, 13), (14, 11)]
qvals8 = [sqrt(i[1] / i[2] - 539/450) / sqrt(8/9) for i in rationals8]
qlist8 = [i[1] for i in rationals8]
plist8 = [i[2] for i in rationals8]
#these are very good surfaces
surfs8 = construct_surfaces(plist8, qlist8, qvals8, prob);
#improvment around tae, but looks ot be worse at r=0.55 ish. and perhaps at r=0.25 ish
plot_surfs(surfs8)
plot_surfs(vcat(surfs, surfs2, surfs3, surfs4, surfs5, surfs6, surfs7, surfs8))
#%%
#add more in the chaotic region, but away from the sepratrix.
rationals9 = [(17, 12), (18, 13)]
qvals9 = [sqrt(i[1] / i[2] - 539/450) / sqrt(8/9) for i in rationals9]
qlist9 = [i[1] for i in rationals9]
plist9 = [i[2] for i in rationals9]
#these are probably the perfect example of more surfaces increasing the jacobian and cooking the spectrum
surfs9 = construct_surfaces(plist9, qlist9, qvals9, prob);
#improvment around tae, but looks ot be worse at r=0.55 ish. and perhaps at r=0.25 ish
plot_surfs(surfs9)
plot_surfs(vcat(surfs, surfs2, surfs3, surfs4, surfs5, surfs6, surfs7, surfs8, surfs9))
#%%
#add surfaces for r >0.6
rationals10 = [(12, 7), (13, 7)]
qvals10 = [sqrt(i[1] / i[2] - 539/450) / sqrt(8/9) for i in rationals10]
qlist10 = [i[1] for i in rationals10]
plist10 = [i[2] for i in rationals10]
#these are probably the perfect example of more surfaces increasing the jacobian and cooking the spectrum
surfs10 = construct_surfaces(plist10, qlist10, qvals10, prob);
#improvment around tae, but looks ot be worse at r=0.55 ish. and perhaps at r=0.25 ish
plot_surfs(surfs10)
plot_surfs(vcat(surfs, surfs2, surfs3, surfs4, surfs5, surfs6, surfs7, surfs8, surfs9, surfs10))
#%%
curr_surfs = vcat(surfs, surfs2, surfs3);
curr_surfs = vcat(surfs, surfs2, surfs3, surfs4);
curr_surfs = vcat(surfs, surfs2, surfs3, surfs4, surfs5);
curr_surfs = vcat(surfs, surfs2, surfs3, surfs4, surfs5, surfs6);
curr_surfs = vcat(surfs, surfs2, surfs3, surfs4, surfs6);
curr_surfs = vcat(surfs, surfs2, surfs3, surfs4, surfs6, [surfs5[2]]);
curr_surfs = vcat(surfs, surfs2, [surfs3[1]], surfs4, surfs6, [surfs5[2]]);
curr_surfs = vcat(surfs, surfs2, surfs3, surfs4, surfs5, surfs6, surfs7);
curr_surfs = vcat(surfs, surfs2, surfs3, surfs4, surfs5, surfs6, surfs7, surfs8, surfs9);
curr_surfs = vcat(surfs, surfs2, surfs3, surfs4, [surfs5[2]], surfs6, surfs7, [surfs8[2]]);
#this batch actually has a tae, think we need to stay away from the sepratrix
curr_surfs = vcat(surfs, surfs2, surfs3, surfs4, surfs7, [surfs8[2]]);
#this is probably the best combination, staying far away from the sepratrix.
curr_surfs = vcat(surfs, surfs2, surfs4, surfs7, [surfs8[2]]);
curr_surfs = vcat(surfs, surfs2, surfs4, surfs7, [surfs8[2]], surfs9);
#think these are the surfaces we want.
curr_surfs = vcat(surfs, surfs2, surfs4, surfs7, [surfs8[2]], surfs10);
total_surfs = vcat(surfs, surfs2, surfs3, surfs4, surfs5, surfs6, surfs7, surfs8, surfs9, surfs10);
length(total_surfs)
plot_surfs(curr_surfs)
plot_surfs(total_surfs)
plot_surfs(surfs)
plot_surfs(vcat(surfs, surfs2, surfs4))
plot_surfs(vcat(surfs, surfs2, surfs4, surfs7))
plot_surfs(vcat(surfs, surfs2, surfs4, surfs7, surfs10))
#%%
surfsl1 = surfs;
surfsl2 = vcat(surfs, surfs2);
surfsl3 = vcat(surfs, surfs2, surfs4);
surfsl4 = vcat(surfs, surfs2, surfs4, [surfs7[1]], [surfs10[1]]);
surfsl5 = vcat(surfs, surfs2, surfs4, surfs7, surfs10, [surfs8[2]]);
#this is hopefully where the `peak' is.
surfsl6 = vcat(surfs, surfs2, surfs3, surfs4, surfs7, surfs10, [surfs8[2]]);
good_rationasl = vcat(rationals, rationals2, rationals3, rationals4, rationals7, [rationals8[2]], rationals10)
surfsl7 = vcat(surfs, surfs2, surfs3, surfs4, surfs5, surfs7, [surfs8[2]], surfs10);
surfsl8 = vcat(surfs, surfs2, surfs3, surfs4, surfs5, surfs6, surfs7, [surfs8[2]], surfs10);
surfsl9 = vcat(surfs, surfs2, surfs3, surfs4, surfs5, surfs6, surfs7, surfs8, surfs10);
surfsl10 = vcat(surfs, surfs2, surfs3, surfs4, surfs5, surfs6, surfs7, surfs8, surfs9, surfs10);

surfs_list = [surfsl1, surfsl2, surfsl3, surfsl4, surfsl5, surfsl6, surfsl7, surfsl8, surfsl9, surfsl10];
#%%
evals, ϕ, ϕft = compute_spectrum_qfm(grids=grids, prob=prob, solver=solver, surfs=curr_surfs);

continuum_plot(evals)

ind = find_ind(evals, 0.2977513)

potential_plot(ϕft, grids, ind)

#%%
#surfsl6 appear to be a good selection of surfaces
length(surfsl6)

save_object("/Users/matt/phd/MID/data/surfaces/chaos_surfs.jld2", surfsl6)
#%%

function compute_flux(prob, grids, surfs)

    #instantiate the grids into arrays. 
    rgrid, θgrid, ζgrid = MID.Structures.inst_grids(grids)

    #initialise the two structs to store the metric and the magnetic field.
    tor_met = MID.Geometry.MetT()
    qfm_met = MID.Geometry.MetT()
    tor_B = MID.Equilibrium.BFieldT()
    qfm_B = MID.Equilibrium.BFieldT()

    #creates the interpolations for the surfaces.
    surf_itp, sd = MID.QFM.create_surf_itp(surfs)

    #compute the gaussian qudrature points for finite elements.
    #ξr, wgr = MID.Construct.FastGaussQuadrature.gausslegendre(grids.r.gp) #same as python!
    #ξθ, wgθ = MID.Construct.FastGaussQuadrature.gausslegendre(grids.θ.gp)
    #ξζ, wgζ = MID.Construct.FastGaussQuadrature.gausslegendre(grids.ζ.gp)

    #struct for storing the intermediate data for the coordinate transform
    CT = MID.QFM.CoordTsfmT()
    #jac = zeros(length(rvals), length(θvals), length(ζvals))
    jac = zeros(length(rgrid), length(θgrid), length(ζgrid))
    djac = zeros(3, length(rgrid), length(θgrid), length(ζgrid))
    B = zeros(3, length(rgrid), length(θgrid), length(ζgrid))

    #for (i, r) in enumerate(rvals), (j, θ) in enumerate(θvals), (k, ζ) in enumerate(ζvals)
    for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)
        MID.QFM.coord_transform!(r, θ, ζ, CT, surf_itp, sd)
        MID.Geometry.toroidal_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], 4.0)
        MID.Equilibrium.compute_B!(tor_B, tor_met, prob.q, prob.isl, prob.isl2, CT.coords[1], CT.coords[2], CT.coords[3])
        MID.QFM.met_transform!(tor_met, qfm_met, CT)
        MID.QFM.B_transform!(tor_B, qfm_B, qfm_met, CT)

        jac[i, j, k] = qfm_met.J[1]
        djac[:, i, j, k] = qfm_met.dJ[:]
        B[:, i, j, k] = qfm_B.B[:]
        #jac_tor[i, j, k] = tor_met.J[1]
        #jac[i, j, k] = qfm_B.B[1]
        #jac_tor[i, j, k] = qfm_B.B[2]
    end
    return B, jac, djac
end
#%%
B, jac, djac = compute_flux(prob, grids, curr_surfs);
#%%
rgrid_plot, θgrid_plot, _ = MID.Structures.inst_grids(grids);
contourf(θgrid_plot, rgrid_plot, jac[:, :, 1], levels=50, title="Jacobian")
contourf(θgrid_plot, rgrid_plot, djac[1, :, :, 1], levels=50, title="dJdr")
contourf(θgrid_plot, rgrid_plot, djac[2, :, :, 1], levels=50, title="dJdt")
#contourf(θgrid_plot, rgrid_plot, B[1, :, :, 1], levels=50, title="B^s")
contourf(θgrid_plot, rgrid_plot, (B[1, :, :, 1]).^2, levels=50, title="B^s^2")
#%%
B2mean = zeros(length(surfs_list));
djacdrmean = zeros(length(surfs_list));
djacdθmean = zeros(length(surfs_list));

for i in 1:length(surfs_list)
    curr_surfs = surfs_list[i]

    B, jac, djac = compute_flux(prob, grids, curr_surfs);

    B2mean[i] = mean(B[1, :, :, :] .^ 2)
    djacdrmean[i] = mean(abs.(djac[1, :, :, :]))
    djacdθmean[i] = mean(abs.(djac[2, :, :, :]))
end
#%%

#ok this actually matches our finding really well

#we may want to change the first set of surfaces, as they artificially look better I think.
#note that Bs is actually ~10 times smaller at 6 vs 1.
#note that Helander et al uses ~60 surfaces, so perhaps we will want more?
nsurfs = [length(surfs) for surfs in surfs_list]
nsurfs[end]
#%%
plot_font = "Computer Modern"
default(fontfamily=plot_font, grid=false, framestyle=:semi, palette=:tol_bright)
cur_colors = palette(:tol_bright);
lfs = 18
xfs = 20
tfs = 16
plot(nsurfs, B2mean, palette=:tol_bright, axis=:left, xlimits=(10, 30), label=L"(B^s)^2", color=cur_colors[1], legendfontsize=lfs, xlabel="Number of QFM Surfaces", xguidefontsize=xfs, xtickfontsize=tfs, ytickfontsize=tfs-3, ylims=(0.0, 0.00005), dpi=1200, legend=:topleft, margin=5Plots.mm, y_foreground_color_border=cur_colors[1], y_foreground_color_axis=cur_colors[1], y_foreground_color_text=cur_colors[1])

#plot!(nsurfs, djacdrmean, axis=:right, palette=:tol_bright, label=L"\partial\mathcal{J}/\partial s")
plot!(twinx(), nsurfs, djacdrmean, palette=:tol_bright, xlimits=(10, 30), axis=:right, linestyle=:dot,linewidth=3, legend=:topright, label=L"|\partial\mathcal{J}/\partial s|", color=cur_colors[2], legendfontsize=lfs, xguidefontsize=xfs, tickfontsize=tfs, dpi=1200, y_foreground_color_border=cur_colors[2], y_foreground_color_axis=cur_colors[2], y_foreground_color_text=cur_colors[2], markerstrokewidth=10)#, yticklabelcolor=cur_colors[2])

plot!(twinx(), nsurfs .+ 100, djacdrmean, xlimits=(10, 30),  palette=:tol_bright, axis=:right, linestyle=:dot,linewidth=1.5, legend=:topright, label=L"|\partial\mathcal{J}/\partial s|", color=cur_colors[2], legendfontsize=lfs, xguidefontsize=xfs, tickfontsize=tfs, dpi=1200, y_foreground_color_border=cur_colors[2], y_foreground_color_axis=cur_colors[2], y_foreground_color_text=cur_colors[2], markerstrokewidth=10)#, yticklabelcolor=cur_colors[2])

savefig("/Users/matt/Dropbox/Apps/Overleaf/Lit Review/pics/QFMPaper/number_of_surfaces.png")
#%%
#plot(nsurfs, djacdθmean)
save_object("/Users/matt/phd/QFMPaper/nsurfs.jld2", nsurfs)
save_object("/Users/matt/phd/QFMPaper/B2mean.jld2", B2mean)
save_object("/Users/matt/phd/QFMPaper/djacdrmean.jld2", djacdrmean)
save_object("/Users/matt/phd/QFMPaper/rationals.jld2", good_rationasl)


