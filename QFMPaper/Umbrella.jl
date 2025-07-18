#copy of the fractal umbrella for our case.
#here we try and actually do Stuarts method to avoid the low order rationals.

using MIDCantori
using MID
using Plots
using LaTeXStrings
#%%
ir1 = [0, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
#either this is somewhat wrong, or the point is actually one of the low order rationals that is wrong...
ir2 = [0, 1, 2, 3, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

MIDCantori.NumberTheory.irrat_from_cf(ir1)
MIDCantori.NumberTheory.irrat_from_cf(ir2)
conv1 = MIDCantori.NumberTheory.convergents(ir1)
#kind of looks like we are going to run out of memory to do the largest cases, guess it will have to be a gadi job?
a1 = conv1[16][2]
b1 = conv1[16][1]

k = 0.0017
@time r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 2, k)
conv2 = MIDCantori.NumberTheory.convergents(ir2)
a2 = conv2[11][2]
b2 = conv2[11][1]
#r02, θ02 = MIDCantori.QFM.anal_periodic_orbit(a2, b2, 4, k)
#%%
#hmm doing this with roots instead has a few that don't converge
#this seems to make the plot looks a bit better.
alist1, blist1 = farey_tree(5, 3, 2, 4, 3)
rats = [(alist1[i], blist1[i]) for i in 1:length(alist1)]


Nvals = length(rats)

kclist = zeros(Nvals);
for i in 1:Nvals
    q = rats[i][1]
    p = rats[i][2]
    kclist[i] = MIDCantori.Residue.find_k_c(q, p)
end
#rats = [(alist[i], blist[i]) for i in 1:length(alist)]
#scatter(blist ./ alist, kclist, markersize=0.9)

kir1 = MIDCantori.Residue.find_k_c(a1, b1)
kir2 = MIDCantori.Residue.find_k_c(a2, b2)

gir1 = surface_guess([(a1, b1)], cantori_q)
gir2 = surface_guess([(a2, b2)], cantori_q)
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
scatter(surface_guess(rats, cantori_q), kclist, markersize=0.9)

scatter!([gir1, gir2], [kir1, kir2])
gl = surface_guess(rats, cantori_q)
ind = 29
scatter!([gl[ind]], [kclist[ind]])
argmin(abs.(gl .-0.605))
rats[29]

MIDCantori.Residue.find_k_c(33, 23)


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

