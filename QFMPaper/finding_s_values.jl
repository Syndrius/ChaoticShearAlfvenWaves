#we need to be able to compare our umbrella, given in terms of rationals or the guess of the rationals
#in terms of the qfm coordinates (s).
#that way we can look at the umbrella and say the rational surface with some (a, b) corresponds to value of s
#then we can more specifically choose our s values for waves to focus on.
using MID
using Plots
using JLD2
#%%

#options
#we could pick an s value and map (s, ϑ=[0, 2π]) to (r, θ) and then modify s until we get the point we want?
#pick the (r, θ) value of the trajectory, i.e. find the periodic orbit, then inverse map to (s, ϑ?) -> probably more objective.
#%%

argmin(abs.(kclist .- 0.0014))
kclist
rats[174]
#rats[191]

a = 115
b = 81
surface_guess([(115, 81)], cantori_q)

r0, θ0 = MIDCantori.QFM.anal_periodic_orbit(a, b, 2, 0.0013)

#%%

surfs = load_object("/Users/matt/phd/MID/data/surfaces/qfm/min_surfs.jld2");
surf_itp, sd = MID.QFM.create_surf_itp(surfs);
CT = MID.QFM.CoordTransformT()

r0[1]


s0, ϑ0, φ0 = MID.Mapping.tor_coords_to_qfm(r0[1], θ0[1], 0.0, CT, surf_itp, sd)

#rats
gl = surface_guess(rats, cantori_q)
scatter(gl, kclist, markersize=0.7)
vline!(surface_guess([(115, 81)], cantori_q))
hline!([0.00115, 0.0013, 0.00145, 0.0017])
vline!([0.558])
vline!([0.61])
surface_guess([(115, 81)], cantori_q)

#so we need to work out what the irrational surfaces are!
b / a
cf = MIDCantori.NumberTheory.continued_fraction(b/a, 10)
conv1 = MIDCantori.NumberTheory.convergents([0, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
conv1[10]
b, a = conv1[10]
vline!(surface_guess([(a, b)], cantori_q))


#%%
alist, blist = farey_tree(5, 3, 2, 4, 3)
#alist, blist = farey_tree(3, 2, 1, 3, 2) #stuart calls this the first 5 levels of the fare tree!
rats = blist ./ alist
perm = sortperm(rats)
rats_s = rats[perm]
alist_s = alist[perm]
blist_s = blist[perm]
γ = 1/2*(1+sqrt(5))

irrats = []
#ok so this process matches Stuart.
for i in 2:length(alist)
    #need to pick which one is more noble...
    ir1 = (blist_s[i] + γ*blist_s[i-1]) / (alist_s[i] + γ*alist_s[i-1]);
    ir2 = (blist_s[i-1] + γ*blist_s[i]) / (alist_s[i-1] + γ*alist_s[i]);
    ir = MIDCantori.NumberTheory.more_noble(ir1, ir2);
    #display(ir)
    push!(irrats, MIDCantori.NumberTheory.continued_fraction(ir, 15))
    #push!(irrats, MIDCantori.NumberTheory.continued_fraction(ir, 15))
end
#display(irrats)
#%%
glir = []
for i in irrats
    conv = MIDCantori.NumberTheory.convergents(i)
    b, a = conv[10]
    push!(glir, surface_guess([(a, b)], cantori_q))
end
#%%
scatter(glir, kclist, markersize=0.7)
vline!(glir[6])

glir[6]
#glir[2] is the right irrat
#glir[5] is the left one
irrats[5]
ir1 = [0, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
ir2 = [0, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

MIDCantori.NumberTheory.irrat_from_cf(ir1)
MIDCantori.NumberTheory.irrat_from_cf(ir2)

#%%
conv1 = MIDCantori.NumberTheory.convergents(ir1)
a1 = conv1[15][2]
b1 = conv1[15][1]

k = 0.0017
r01, θ01 = MIDCantori.QFM.anal_periodic_orbit(a1, b1, 4, k)
conv2 = MIDCantori.NumberTheory.convergents(ir2)
a2 = conv2[12][2]
b2 = conv2[12][1]
r02, θ02 = MIDCantori.QFM.anal_periodic_orbit(a2, b2, 4, k)
#%%

N = 2

isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isls = MID.IslandT[isl1, isl2]
geo = init_geo(R0=1.0)
prob = init_problem(q=cantori_q, geo=geo, met=:cylinder, isls=isls, type=:flux)
#%%

#cool, two flux surfaces to focus on!
poincare_plot(prob, 1000, [r01[1], r02[1]], [θ01[1], θ02[1]], ylimits=(0.5, 0.6667))
#%%
s01, ϑ01, _ = MID.Mapping.tor_coords_to_qfm(r01[1], θ01[1], 0.0, CT, surf_itp, sd)
s02, ϑ02, _ = MID.Mapping.tor_coords_to_qfm(r02[1], θ02[1], 0.0, CT, surf_itp, sd)
#%%
#now we need to be able to create an s grid that is as close to our two values as possible.
#k1 -> s1=0.055806696, s2=0.6115465
#k2 -> s1=0.55805786, s2=0.61007268
#k3 -> s1=0.55810525, s2=0.60871189
#k4 -> s1=0.55808634, s2=0.606684930268
#so it changes, guess we can just stick with our original grid.
