#need to figure out if the different surfaces result in a different mapping in each case
#maybe we can justify our weird choice of flux surface?
#this shows our results are even more cooked than expected.
using MID
using JLD2
using MIDCantori
using Plots; plotlyjs()
using MIDViz
#%%

#ooft, so there is a significant difference in the third decimal place
#this is enough to shift us a flux surface or even a few.
#think we need to actually get the irrational surface...
#neeed to fix the fkn umbrella.

alist, blist = farey_tree(6, 4, 3, 3, 2)
rats = [(alist[i], blist[i]) for i in 1:length(alist)]
kclist = MIDCantori.Residue.umbrella(rats)
gl = surface_guess(rats, cantori_q)
scatter(gl, kclist, title="k0")
#%%
fl, _ = cantori_q(0.5417)
cf = MIDCantori.NumberTheory.continued_fraction(fl, 15)
conv = MIDCantori.NumberTheory.convergents(cf)

#(48, 35) will work as second example, for k05 it is v good 
#k08 there is a decent example, but we can see the onset of chaos
#k10
surface_guess([(48, 35)], cantori_q)
#%%

surfs = load_object("/Users/matt/phd/MID/data/surfaces/qfm/k11_surfs.jld2");
surf_itp, sd = MID.create_surf_itp(surfs);
ind1 = find_ind(gl, 0.607)
ind2 = find_ind(gl, 0.5416)
#MIDCantori.Residue.find_k_c(62, 45)
rp1, θp1 = periodic_orbit(rats[ind1]..., 8, 0.0011)
#rp2, θp2 = periodic_orbit(rats[ind2]..., 8, 0.0013)
rp2, θp2 = periodic_orbit(48, 35, 8, 0.0011)
#%%
r = rp1[1]
#r = 0.59259
θ = θp1[1]
ζ = 0.0

MID.Mapping.tor_coords_to_qfm(r, θ, ζ, CT, surf_itp, sd)
r = rp2[1]
#r = 0.59259
θ = θp2[1]
#θ = 0.0

MID.Mapping.tor_coords_to_qfm(r, θ, ζ, CT, surf_itp, sd)
#%%

CT = MID.QFM.CoordTransformT();

s = 0.54040


MID.QFM.coord_transform!(s, 0.0, 0.0, CT, surf_itp, sd);
CT.coords

#%%
#perhaps we try and map the entire guess list, to see roughly where the regions are cooked af.
#but i think this is just further evidence that this doesn't work.
surfs = load_object("/Users/matt/phd/MID/data/surfaces/qfm/k11_surfs.jld2");
surf_itp, sd = MID.create_surf_itp(surfs);

gl11 = zeros(length(gl))
for i in 1:length(gl)
    sm, ϑm, φm = MID.Mapping.tor_coords_to_qfm(gl[i], 0.0, 0.0, CT, surf_itp, sd)
    gl11[i] = sm
end

scatter(gl11, kclist, title="k11")

#these do not seem to match the guess list...
#this is because the gl values are assuming unperturbed, so it is mapping a sort of random value tomewhere else.

#%%
#so the two methods of transforming the coordinates yeild v different results...
#not ideal.
k=0.0012
geo = init_geo(R0=1.0)
isl1 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isl2 = init_island(m0=4, n0=-3, A=k/4, flux=true)

isls = MID.Geometry.IslandT[isl1, isl2]

prob = init_problem(geo=geo, q=cantori_q, isls=isls, met=:cylinder, type=:flux)

Ntraj = 100;
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
ψlist = collect(LinRange(0.4, 0.8, Ntraj));
Nlaps = 500;

#would be good if we could plot the poincare with like a gradient of colours to see where they end up.
poincare_plot(prob, Nlaps, ψlist, zeros(Ntraj), ylimits=(0.5, 0.67))
#%%
poincare_plot(prob, Nlaps, ψlist, zeros(Ntraj), surfs=surfs, ylimits=(0.5, 0.67))
