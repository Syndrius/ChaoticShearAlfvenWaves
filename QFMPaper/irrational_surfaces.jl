#finding the irrational value of our two surfaces (or just one! since the other is cooked)

using MID
using MIDCantori
using MIDViz
using Plots; gr()
using Plots; plotlyjs()
#%%
alist, blist = farey_tree(6, 4, 3, 3, 2)
rats = [(alist[i], blist[i]) for i in 1:length(alist)]
kclist = MIDCantori.Residue.umbrella(rats)
gl = surface_guess(rats, cantori_q)
scatter(gl, kclist)
#%%
indl = find_ind(gl, 0.608)
indr = find_ind(gl, 0.607)
rats[indl]
rats[indr]
#%%
γ = 1/2 * (1+sqrt(5))

#ir1 = 1.4362215818857254
ir1 = (rats[indl][1]+γ*rats[indr][1]) / (rats[indl][2] + γ*rats[indr][2])
ir2 = (rats[indr][1]+γ*rats[indl][1]) / (rats[indr][2] + γ*rats[indl][2])

cf1 = MIDCantori.NumberTheory.continued_fraction(ir1, 15)
cf2 = MIDCantori.NumberTheory.continued_fraction(ir2, 15)

conv1 = MIDCantori.NumberTheory.convergents(cf1)
conv2 = MIDCantori.NumberTheory.convergents(cf2)

#%%

rat1 = conv1[10]
rat2 = conv2[10]


MIDCantori.Residue.find_k_c(rat1...) #looks like this one is slightly better, not certain it is the best but whatever.
MIDCantori.Residue.find_k_c(rat2...)
#%%
k = 0.0012
geo = init_geo(R0=1.0)
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)

isls = MID.Geometry.IslandT[isl1, isl2]

prob = init_problem(geo=geo, q=cantori_q, isls=isls, met=:cylinder, type=:flux)
#%%
Ntraj = 350;
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
ψlist = collect(LinRange(0.5, 0.7, Ntraj));
Nlaps = 500;

#would be good if we could plot the poincare with like a gradient of colours to see where they end up.
x, z = poincare_plot(prob, Nlaps, ψlist, zeros(Ntraj), ylimits=(0.5, 0.67))
#%%
rp, θp = MIDCantori.Orbits.anal_periodic_orbit(rat1..., 8, k)
r1, θ1 = poincare_plot(prob, Nlaps, [rp[1]], [θp[1]])
#%%
#good enough, looks like the correct flux surface, however, we may need to change this!
#think this actually works pretty well, flux surface is clearly a cantori for k=0.0015, but pretty heked for k=0.0017
#so this should give us a pretty good idea.
poincare_plot(x, z, ylimits=(0.575, 0.63), color=:black)
poincare_plot(r1, θ1, overlay=true, markersize=2)


#%%
#now lets consider e
