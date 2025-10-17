
#creation of the umbrella pic for aapps.
#main difference between paper is only a single flux surface is highlighted.
using MIDCantori
using JLD2
using Plots; gr()
#using Plots; plotlyjs()
using MID
using LaTeXStrings
#%%
ylims = (0.0, 1.5e-3)
#%%

#clearly got some nonsense in it.
kc_fir109 = load_object("/Users/matt/phd/MID/QFMPaper/results/umbrella/umb_farey_ir109.jld2")
#also got some garbage.
kc_fir88 = load_object("/Users/matt/phd/MID/QFMPaper/results/umbrella/umb_farey_ir88.jld2")

kc_f244 = load_object("/Users/matt/phd/MID/QFMPaper/results/umbrella/umb_farey244.jld2")

kc_f255 = load_object("/Users/matt/phd/MID/QFMPaper/results/umbrella/umb_farey255.jld2")


kc_ir_rat47 = load_object("/Users/matt/phd/MID/QFMPaper/results/umbrella/umb_ir_rat47.jld2")
#%%

ratsfir109 = umbrella_rats((4, 3), (3, 2), 10, 9, :farey_ir)
ratsfir88 = umbrella_rats((4, 3), (3, 2), 8, 8, :farey_ir)
ratsf244 = umbrella_rats((4, 3), (3, 2), 4, 4, :farey2)
ratsf255 = umbrella_rats((4, 3), (3, 2), 5, 5, :farey2)
ratsir_rat = umbrella_rats((4, 3), (3, 2), 4, 7, :ir_rat)
dvsfir109 = [i[1] / i[2] for i in ratsfir109]
dvsfir88 = [i[1] / i[2] for i in ratsfir88]
dvsf244 = [i[1] / i[2] for i in ratsf244]
dvsf255 = [i[1] / i[2] for i in ratsf255]
dvsir_rat = [i[1] / i[2] for i in ratsir_rat]
#%%
scatter(dvsfir109, kc_fir109)
scatter(dvsfir88, kc_fir88, ylimits=ylims)
scatter(dvsf244, kc_f244, ylimits=ylims)
scatter(dvsf255, kc_f255, ylimits=ylims)
scatter(dvsir_rat, kc_ir_rat47, ylimits=ylims)
scatter!(dvsf244, kc_f244, ylimits=ylims)
#%%
#kind of the same as before
#think farey 44 is best, just needs some pruning.

scatter(dvsf244, kc_f244, ylimits=ylims)

find_ind(dvsf244, 1.424658)
to_remove = [236, 87, 86, 88, 67, 82, 74]
#%%
kcs = []
dvs = []
for i in 1:length(kc_f244)
    if i in to_remove
        continue
    else
        push!(kcs, kc_f244[i])
        push!(dvs, dvsf244[i])
    end
end
#%%
scatter(dvs, kcs)
#%%
plot_font = "Computer Modern"
#default(fontfamily=plot_font, grid=false, framestyle=:semi, palette=:tol_bright)
default(fontfamily=plot_font, grid=false, palette=:tol_bright)
cur_colors = palette(:tol_bright);
#lfs = 16
xfs = 28
yfs = 26
ytfs = 10
xtfs = 14
lfs = 16
xlab = L"\psi"
ylab = L"k_c"
conv_num = 15 #probably want this to be 15 for actual results, finding the all the orbits takes ~20 mins, smaller number for changing the plot vars.
ylms = (0.0005, 0.0015)
xlms = (4/3, 3/2)
tfs = 16
dpi = 1200
msize = 0.7
msize1 = 1.7
msize2 = 1.7
#1 is used for the (6, 4) branch!
c1 = cur_colors[2]
c2 = cur_colors[3]
alval = 0.5
#c1 =:black
#c2 =:black

shape1 = :o
shape2 = :o
xtks = surface_guess([(3, 2), (7, 5), (4, 3)], cantori_q)
xtks = [3/2, 7/5, 4/3, 10/7, 11/8]
#xtks = [xtks ; gir1; gir2]
#xtkslab = [L"3/2", L"7/5", L"4/3", L"\iota_1", L"\iota_2"]
xtkslab = [L"3/2", L"7/5", L"4/3", L"10/7", L"11/8"]
#%%
conv_num = 7
#may not be perf, especially since it is not very noble!
#but it is in the correct spot! so its fine!
ir1 = [1, 2, 1, 2, 2, 1, 5, 1, 5, 1, 5, 1, 5, 1, 5, 1, 5, 1, 5, 1, 5, 1, 5, 1] #not very noble!
ir2 = [1, 2, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
conv1 = MIDCantori.NumberTheory.convergents(ir1)
conv2 = MIDCantori.NumberTheory.convergents(ir2)
a1 = conv1[conv_num][1]
b1 = conv1[conv_num][2]
a2 = conv2[conv_num][1]
b2 = conv2[conv_num][2]

kcir1 = find_k_c(a1, b1)
kcir2 = find_k_c(a2, b2)

#%%
scatter(dvs, kcs, markersize=1.0, label=false, xlimits=xlms, ylimits=ylms, dpi=dpi, ylabel=ylab, xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, legendfontsize=lfs, xticks=(xtks, xtkslab), color=:black, left_margin=3*Plots.mm, right_margin=3Plots.mm)
scatter!([a1/b1], [kcir1], color=c1, markershape=:dtriangle, markersize=5.0, label=false)
scatter!([a2/b2], [kcir2], color=c2, markershape=:dtriangle, markersize=5.0, label=false)
#hline!([1.3e-3], label=L"k=0.0013", color=:grey, alpha=0.6, linewidth=2)

savefig("/Users/matt/phd/MID/QFMPaper/results/umbrella.png")

