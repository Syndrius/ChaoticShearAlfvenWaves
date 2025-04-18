
#here we generate neat copies of the poincare plots and surfaces needed for puplication.
using MID
using MIDviz
#%%
#save_dir = "/Users/matt/Dropbox/Apps/Overleaf/Lit Review/pics/QFMPaper/benchmark/"
save_dir = "/Users/matt/Dropbox/Apps/Overleaf/Lit Review/pics/QFMPaper/chaos/"

#%%
R0=4.0

#amp needs further thought!
#define the non-resonant island
#with chaos_q, k=0.0025 is very chaotic while sitll having an inner island chain
#0.0027 seems to be about when there is no structures left
#think k=0.0022 may be the best bet, there are very clear structures but it is also very chaotic elsewhere.
k = 0.0006
isl = init_island(m0=3, n0=-2, A=k/3)
isl2 = init_island(m0=4, n0=-3, A=k/4)
#k = 0.008
#isl = init_island(m0=3, n0=2, A=k/3)
#isl2 = init_island(m0=4, n0=-3, A=0.0)

geo = init_geo(R0=R0)

#to solve non-Hermitian
#flr = MID.Structures.FLRT(δ = 1e-18)

prob = init_problem(q=qfm_q, geo=geo, isl=isl, isl2=isl2)#, flr=flr)
#%%
#Define parameters needed for the poincare plot
Ntraj = 200;
rlist = collect(LinRange(0.3, 0.8, Ntraj));
#rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

#poincare_plot(qfm_benchmark_q, slab_to_plot, Nlaps, Ntraj, 0, 0.0, 0.0, 0.0, R0, isl, MID.Structures.no_isl, rlist)
x, z = poincare_plot(prob, Nlaps, Ntraj, rlist,  R0)#, filename=save_dir * "original_poincare.png");

#%%
#bit unsure about the limits for this one tbh.
MIDViz.plot_poincare(x, z, axis=true, yguidefontsize=25, xguidefontsize=20, palette=:tol_bright, ylimits=(0.4, 0.7), filename=save_dir*"poincare.png")

#%%

qlist, plist = farey_tree(4, 1, 1, 2, 1)
#wrong q profile for this now!
guess_list = sqrt(0.5) .* sqrt.(qlist ./ plist .- 0.99);
@time surfs = construct_surfaces(plist, qlist, guess_list, prob);
#%%
plot_surfs(surfs, xguidefontsize=20, yguidefontsize=25, legend=false, palette=:tol_bright, filename=save_dir*"surfaces.png")
#%%

Ntraj = 40;
#rlist = collect(LinRange(0.4, 0.65, Ntraj));
#so r < 0.4 seems to be completly cooked for every situation...
#not ideal
rlist = collect(LinRange(0.4, 0.8, Ntraj));
Nlaps = 500;

#so this doesn't work v nice.
#guess this tells us that the surfaces are cooked

#this looks to sort of be working now, probably need many more surfaces tbh.
#poincare_plot(qfm_benchmark_q, slab_to_plot, Nlaps, Ntraj, 0, 0.0, 0.0, 0.0, R0, isl, MID.Structures.no_isl, rlist)
xq, zq = poincare_plot(prob, Nlaps, Ntraj, rlist,  R0, surfs=surfs)#, filename=save_dir * "qfm_poincare.png");
#%%
#bit unsure about the limits for this one tbh.
MIDViz.plot_poincare(xq, zq, axis=true, yguidefontsize=25, xguidefontsize=20, ylimits=(0.4, 0.7), paletta=:tol_bright, ylabel=L"s", xlabel=L"\vartheta", filename=save_dir*"qfm_poincare.png")

#%%

