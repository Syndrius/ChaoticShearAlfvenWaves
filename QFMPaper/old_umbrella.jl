#copy of the fractal umbrella for our case.
#here we try and actually do Stuarts method to avoid the low order rationals.
#increasing the depth just created a case where it gets cooked, 
#either need to be happy enough with the curretn state and just fix it up a bit
#or we need to improve the solver/integrator.

using MIDCantori
using MID
using Plots
using LaTeXStrings
#%%
ir1 = [0, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
#either this is somewhat wrong, or the point is actually one of the low order rationals that is wrong...
#looks like this is a bit off!
#may have fixed this now
#but now we have changed to the most resilaint flux surface at ~0.0592 ish.
ir2 = [0, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

ir1_val = MIDCantori.NumberTheory.irrat_from_cf(ir1)
ir2_val = MIDCantori.NumberTheory.irrat_from_cf(ir2)
conv1 = MIDCantori.NumberTheory.convergents(ir1)
#kind of looks like we are going to run out of memory to do the largest cases, guess it will have to be a gadi job?
a1 = conv1[16][2]
b1 = conv1[16][1]

k = 0.0017
#@time r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 2, k)
conv2 = MIDCantori.NumberTheory.convergents(ir2)
a2 = conv2[16][2]
b2 = conv2[16][1]
#r02, θ02 = MIDCantori.QFM.anal_periodic_orbit(a2, b2, 4, k)
#%%
res_file = "results/test.txt"
f = open(res_file)

kclist_str = readline(f)
kc_split = split(replace(kclist_str[2:end-1], ","=> ""))
kclist = parse.(Float64, kc_split)
close(f)
#%%
#rat1 = (4, 3)
#rat2 = (3, 2)
#depth = 5
#new_rats = MIDCantori.Residue.create_umbrella_rats(rat1, rat2, depth, 2:10, 10)
#this method works pretty well
#still afew odd points, but nothing serious.
#somehow this only took ~5 minutes (with 48 cores!)
#meaning we can probably run this again with a deaper farey tree.
alist, blist = farey_tree(4, 3, 2, 4, 3)

rats = [(alist[i], blist[i]) for i in 1:length(alist)]

dvs = [i[2]/i[1] for i in rats]

perm = sortperm(dvs)

#only need to sort the global one I think!
rats_sort = rats[perm]
g_rats = Tuple{Int64, Int64}[]
for i in 1:length(rats_sort) - 1
    r1 = rats_sort[i]
    r2 = rats_sort[i+1]
    al, bl = farey_tree(4, r1[1], r1[2], r2[1], r2[2])
    rl = [(al[j], bl[j]) for j in 3:length(al)]
    g_rats = [g_rats ; rl]
end

gl = surface_guess(g_rats, cantori_q)
#%%
#first irrational!
argmax(kclist[140:end])
gl[140+6]
kclist[140:end][6]

scatter(gl, kclist, markersize=0.9)
vline!([gl[140+6]])

#second inrrat
new_kclist = kclist
new_kclist[74] = 0.0
argmax(new_kclist)

kclist[84]
g_rats[84]
gl[84]
scatter(gl, kclist, markersize=0.9)
vline!([gl[84]])
hline!([0.00145], label="k3")
hline!([0.00115], label="k1")
hline!([0.0013], label="k2")
hline!([0.00125], label="k15")
#%%

#doesn't work atm, shows that one of our irrationals is wrong!
#kir1 = MIDCantori.Residue.find_k_c(a1, b1)
#kir2 = MIDCantori.Residue.find_k_c(a2, b2)

gir1 = surface_guess([(a1, b1)], cantori_q)
#irrat 2 is not working at all atm!
#have changed this now!
gir2 = surface_guess([(a2, b2)], cantori_q)
#don't think this is neede now, new one is 0.610427
#this value is not really correct but it is kinda close enough
gir2 = 0.5915
#hline!([0.0014])
#%%
#now we actually plot the damn thing!
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
ylms = (0.0, 0.0018)
xlms = (0.5, 2/3)
tfs = 16
dpi = 1200
msize = 0.7
msize1 = 1.7
msize2 = 1.7
c1 = cur_colors[1]
c2 = cur_colors[2]
alval = 0.5
#c1 =:black
#c2 =:black

shape1 = :o
shape2 = :o

xtks = surface_guess([(3, 2), (7, 5), (4, 3)], cantori_q)
xtks = [xtks ; gir1; gir2]
xtkslab = [L"3/2", L"7/5", L"4/3", L"\iota_1", L"\iota_2"]
#%%
#this is going to be a busy plot.
scatter(gl, kclist, markersize=0.9, label=false, ylimits=ylms, xlimits=xlms, dpi=dpi, ylabel=ylab, xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=xtfs, ytickfontsize=ytfs, legendfontsize=lfs, xticks=(xtks, xtkslab))

savefig("/Users/matt/phd/MID/QFMPaper/results/base_umbrella.png")
vline!([gir1], color=cur_colors[1], alpha=alval, linestyle=:dash, label=false)
vline!([gir2], color=cur_colors[2], alpha=alval, linestyle=:dash, label=false)
#vline!([gir1], label=L"\iota_1", color=cur_colors[1], alpha=0.3, linestyle=:dash)
#vline!([gir2], label=L"\iota_2", color=cur_colors[2], alpha=0.3, linestyle=:dash)

savefig("/Users/matt/phd/MID/QFMPaper/results/irr_umbrella.png")

#unsure how to label this!
#we may want the irratioanls to be points?
#these could be axis labels like Stuart, same as for the irrats.
#these lines show us that we need more data points in the umbrella!
hline!([0.0005], color=:grey, alpha=alval, label=false)
hline!([0.0011], color=:grey, alpha=alval, label=false)
hline!([0.0013], color=:grey, alpha=alval, label=false)
hline!([0.0015], color=:grey, alpha=alval, label=false)
hline!([0.0017], color=:grey, alpha=alval, label=false)

savefig("/Users/matt/phd/MID/QFMPaper/results/umbrella.png")

#looks like we will need more points towards our key rationals!


#scatter!([gir1, gir2], [kir1, kir2])
#gl = surface_guess(rats, cantori_q)
#ind = 29
#scatter!([gl[ind]], [kclist[ind]])
#argmin(abs.(gl .-0.605))
#rats[29]

#MIDCantori.Residue.find_k_c(33, 23)
#ir_gl = surface_guess([(a1, b1), (a2, b2)], cantori_q)
#vline!(ir_gl)
#%%
#think we need to create our umbrella more like
#consider a farey tree, then between each rational pair, take the most noble irrational,
#i.e. a_1+γ a_2) / (b_1+γ b_2) or whatever
#then consider the convergents of that irrational, hopefully these are not to large, however, we may still need to more accuratly compute the residue, by changing the fidl line tracing or the nlsolve.
#but perhaps this is something we need to talk about with Zhisong.
alist, blist = farey_tree(7, 3, 2, 4, 3)

rats = [(alist[i], blist[i]) for i in 1:length(alist)]

dvs = [i[2]/i[1] for i in rats]

perm = sortperm(dvs)

#only need to sort the global one I think!
rats_sort = rats[perm]
g_rats = Tuple{Int64, Int64}[]
γ = 1/2*(1+sqrt(5))
for i in 1:length(rats_sort) - 1
    rat1 = rats_sort[i]
    rat2 = rats_sort[i+1]
    ir1 = (rat1[2] + rat2[2]*γ) / (rat1[1] + rat2[1]*γ)
    ir2 = (rat2[2] + rat1[2]*γ) / (rat2[1] + rat1[1]*γ)
    cf1 = MIDCantori.NumberTheory.continued_fraction(ir1, 15)
    cf2 = MIDCantori.NumberTheory.continued_fraction(ir2, 15)
    convs1 = MIDCantori.NumberTheory.convergents(cf1)
    convs2 = MIDCantori.NumberTheory.convergents(cf2)
    display(convs1[10])
    display(convs2[10])
    a1 = convs1[10][2]
    b1 = convs1[10][1]
    a2 = convs2[10][2]
    b2 = convs2[10][1]
    push!(g_rats, (a1, b1))
    push!(g_rats, (a2, b2))
end

length(g_rats)

unique(g_rats)

#%%
#again, it sort of looks alright, hard to tell if this approach is better or worse tbh!
#this is sort of just like the other appraoches, good bits and bad bits
#this was very quick at least.
#but there is clearly some nonsense
#which means we probably need more accurate field line tracing!
#and or numerical solving.
#both are shithouse.
gl = surface_guess(g_rats, cantori_q)

scatter(gl, kclist, markersize=0.7, ylimits=(0.0, 0.0018))


############################################################
#%%
#now testing removing them later.
#much better!
#so some of these are still just wrong....
#think going up to 8 is probably fine.
#maybe 9 for final copy?
alist2, blist2 = farey_tree(4, 3, 2, 4, 3) #used to remove the lowest order rationals.
rats2 = [(alist2[i], blist2[i]) for i in 1:length(alist2)]

rats_p = Tuple{Int64, Int64}[]
kc_p = []
for i in 1:length(rats)
    if !(rats[i] in rats2)
        push!(rats_p, rats[i])
        push!(kc_p, kclist[i])
    end
end

p = scatter(surface_guess(rats_p, cantori_q), kc_p, markersize=0.9)

#%%

hl_p = []
for i in rats2
    push!()
end

gl_p = surface_guess(rats2, cantori_q)

for i in gl_p
    vline!([i])
end
display(p)
#%%
#test results from gadi!
#looks even worse not ideal.
#think we probably need to see if we can get the right behaviour for a single rational
#then expand it.
f = open("results/test.txt")
f = open("test.txt")

kclist_str = readline(f)
kc_split = split(replace(kclist_str[2:end-1], ","=> ""))
kclist = parse.(Float64, kc_split)
close(f)
#%%
#rat1 = (4, 3)
#rat2 = (3, 2)
#depth = 5
#new_rats = MIDCantori.Residue.create_umbrella_rats(rat1, rat2, depth, 2:10, 10)
#this method works pretty well
#still afew odd points, but nothing serious.
#somehow this only took ~5 minutes (with 48 cores!)
#meaning we can probably run this again with a deaper farey tree.
alist, blist = farey_tree(4, 3, 2, 4, 3)

rats = [(alist[i], blist[i]) for i in 1:length(alist)]

dvs = [i[2]/i[1] for i in rats]

perm = sortperm(dvs)

#only need to sort the global one I think!
rats_sort = rats[perm]
g_rats = Tuple{Int64, Int64}[]
for i in 1:length(rats_sort) - 1
    r1 = rats_sort[i]
    r2 = rats_sort[i+1]
    al, bl = farey_tree(4, r1[1], r1[2], r2[1], r2[2])
    rl = [(al[j], bl[j]) for j in 3:length(al)]
    g_rats = [g_rats ; rl]
end

gl = surface_guess(g_rats, cantori_q)

scatter(gl, kclist, markersize=0.7, ylimits=(0.0, 0.0015))
#%%
#new aproach to umbrella.
#find the ~17 lowest rationals, then do a farey tree between each one!
#looks like this could work, but will take ages, naturally.
#so ~270 ish took 14 hours on 48 procs.
#try to get a similar number?
alist, blist = farey_tree(5, 3, 2, 4, 3)
length(alist)

rats = [(alist[i], blist[i]) for i in 1:length(alist)]

dvs = [i[2]/i[1] for i in rats]

perm = sortperm(dvs)

#only need to sort the global one I think!
rats_sort = rats[perm]
g_rats = Tuple{Int64, Int64}[]
for i in 1:length(rats_sort) - 1
    r1 = rats_sort[i]
    r2 = rats_sort[i+1]
    al, bl = farey_tree(4, r1[1], r1[2], r2[1], r2[2])
    rl = [(al[j], bl[j]) for j in 3:length(al)]
    g_rats = [g_rats ; rl]
end

length(g_rats)

scatter(surface_guess(g_rats, cantori_q))

g_rats

#al, bl = farey_tree(4, 24, 17, 25, 18)

#rl = [(al[i], bl[i]) for i in 4:length(al)]
#%%
kcl = zeros(length(g_rats))
for i in 1:length(g_rats)
    kcl[i] = MIDCantori.Residue.find_k_c(g_rats[i][1], g_rats[i][2])
end
#%%
gl = surface_guess(g_rats, cantori_q)
scatter(gl, kcl, markersize=0.7)

surface_guess([(7, 5)], cantori_q)
