
#see if we can compute the flux ∫ B^s Jdϑdφ, and see how it changes as the number of surfaces changes
#the flux will decrease with more surfaces
#but we expect the Jacobian to become more volatile with more surfaces,
#ideally, there will be a sweet spot, where the flux is getting diminishing returns but the jacobian is yet to go bananas.
#looks like we can just compute B^s and take the average for a fixed s. 
#should approx the flux, probably worthwhile plotting as a function of s,
#expect we will see the flux will be further from zero near sepratrix etc.
#hopefully we can see behaviour we want.
#once we have determined a good set of surfaces, we probably should check if changing the M, N, MM params changes the surfaces much.

using MID
using MIDViz
using Plots; plotlyjs()
using JLD2
using Statistics
#%%
R0=4.0

geo = init_geo(R0=R0)
k = 0.0006
isl = init_island(m0=3, n0=-2, A=k/3)
isl2 = init_island(m0=4, n0=-3, A=k/4)

prob = init_problem(q=qfm_q, geo=geo, isl=isl, isl2=isl2)#, flr=flr)
#%%
qlist3, plist3 = farey_tree(3, 1, 1, 5, 2)
qlist4, plist4 = farey_tree(4, 1, 1, 5, 2)
qlist5, plist5 = farey_tree(5, 1, 1, 5, 2)
#guess_list = 0.5 .* ones(length(qlist));
guess_list3 = @. sqrt(qlist3 / plist3 - 0.95) / sqrt(1.5);
guess_list4 = @. sqrt(qlist4 / plist4 - 0.95) / sqrt(1.5);
guess_list5 = @. sqrt(qlist5 / plist5 - 0.95) / sqrt(1.5);
#needed for the single island chain case, this one is probably at the sepratrix I guess.
@time surfs3 = construct_surfaces(plist3, qlist3, guess_list3, prob);
@time surfs4 = construct_surfaces(plist4, qlist4, guess_list4, prob);
#surfs 5 is essentially useless, takes ages, and mainly adds surfaces for large r, where there arealready enough
@time surfs5 = construct_surfaces(plist5, qlist5, guess_list5, prob);
plot_surfs(surfs3)
#%%
extra_surfs = construct_surfaces([15, 13, 20, 12], [16, 15, 21, 29], [0.3, 0.3, 0.26, 0.95], prob);

#%%
#the farey tree is not giving any in the chaotic region lol
@time surfsc = construct_surfaces([8, 7], [11, 10], [0.5, 0.5], prob);
#still no surfaces near the island at r=0.45!
#the numbers in here clearly show that the farey tree is a stupid way of doing this
#choosing all combinations of lowest rationals would pick these up, so we need a combination of the two.
#perhaps we should just say we choose surfaces with the lowest rationals, and get extras using mediats
#and randomly pick extras to fit in gaps.
@time surfsc2 = construct_surfaces([3, 4, 10], [4, 5, 11], [0.5, 0.5, 0.25], prob);
@time surfsc3 = construct_surfaces([5, 30], [6, 31], [0.35, 0.1], prob);
#(29, 20) is another valid option, but seems to be unable to find.
#(22, 15) is another valid option, but seems to be unable to find.
#(13, 9) did not actually take very long!
@time surfsc4 = construct_surfaces([10, 9], [13, 13], [0.42, 0.52], prob);
#for some reason we didin't already have (7, 5)! these are meant to fill in small gaps in chaotic region!
#also add (16, 13)
#and (14, 11)
@time surfsc5 = construct_surfaces([5, 13, 11], [7, 16, 14], [0.492, 0.37, 0.38], prob);
plot_surfs(surfsc2)
plot_surfs(surfsc3)
plot_surfs(surfsc4)
plot_surfs(surfsc5)
#%%
save_object("surfs3.jld2", surfs3)
save_object("surfs4.jld2", surfs4)
save_object("surfs5.jld2", surfs5)
save_object("extra_surfs.jld2", extra_surfs)
#arr = [1/1, 26/25, 25/23, 11/10, 7/6, 6/5, 5/4, 4/3, 11/8, 7/5, 10/7, 3/2, 8/5, 5/3, 7/4, 9/5, 11/6, 2/1, 13/6, 11/5, 9/4, 7/3, 12/5, 5/2]
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
    ξr, wgr = MID.Construct.FastGaussQuadrature.gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = MID.Construct.FastGaussQuadrature.gausslegendre(grids.θ.gp)
    ξζ, wgζ = MID.Construct.FastGaussQuadrature.gausslegendre(grids.ζ.gp)

    #struct for storing the intermediate data for the coordinate transform
    CT = MID.QFM.CoordTsfmT()
    #jac = zeros(length(rvals), length(θvals), length(ζvals))
    jac = zeros(grids.r.N, grids.θ.N, grids.ζ.N)
    djac = zeros(3, grids.r.N, grids.θ.N, grids.ζ.N)
    B = zeros(3, grids.r.N, grids.θ.N, grids.ζ.N)

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
Nr = 50
Nθ = 30
Nζ = 9
rgrid = init_grid(type=:rf, N = Nr, start = 0.05, stop =0.99, sep1=0.4, sep2=0.6, frac=0.5)
θgrid = init_grid(type=:af, N = Nθ, pf=3)
ζgrid = init_grid(type=:af, N = Nζ, pf=-2)

rgrid_plot, θgrid_plot, _ = MID.Structures.inst_grids(grids);
grids = init_grids(rgrid, θgrid, ζgrid)
#%%
#going from 4 to 3 (with all the others!) is a huge improvment
#as b^s near r=0.55 is srastically reduced due to the extra surfaces there
#the Jacobian doesn't get worse, and the extreme parts of the jacobain are more contained to the chaotic region.
#with 5, the dJdr goes to ~750 and Bfield goes to 0.0075
#this is worse on all fronts, so 5 is not good, it may perhaps be the case that some of 5 is good though.
qarr = @. sqrt(qlist5 / plist5 - 239/240) / sqrt(5/3)
#written as (q, p) pairs
defos = [(1, 1), (5, 2), (4, 3), (3, 2), (5, 4)]
defqs = [sqrt(pair[1] / pair[2] - 239/240) / sqrt(5/3) for pair in defos]
defqlist = [pair[1] for pair in defos]
defplist = [pair[2] for pair in defos]
#def_surfs = 
#using all of the surfaces in 4 may not be required.
#this is a very basic spread of the surfaces, that should be our starting point, if we were to do this in a better way.
plist_extra = [15, 13, 20, 12]
qlist_extra = [16, 15, 21, 29]
plist_surfsc = [8, 7]
qlist_surfsc = [11, 10]
plist_surfsc2 = [3, 4, 10]
qlist_surfsc2 = [4, 5, 11]
plist_surfsc3 = [5, 30]
qlist_surfsc3 = [6, 31]
plist_surfsc4 = [10, 9]
qlist_surfsc4 = [13, 13]
plist_surfsc5 = [5, 13, 11]
qlist_surfsc5 = [7, 16, 14]
(qlist5[10], plist5[10])
total_q = vcat(qlist5, qlist_extra, qlist_surfsc, qlist_surfsc2, qlist_surfsc3, qlist_surfsc4, qlist_surfsc5);
total_p = vcat(plist5, plist_extra, plist_surfsc, plist_surfsc2, plist_surfsc3, plist_surfsc4, plist_surfsc5);
total_surfs = vcat(surfs5, extra_surfs, surfsc, surfsc2, surfsc3, surfsc4, surfsc5);
ind = 39 #33 and 39 is a double up!
(total_q[ind], total_p[ind])
#(q, p) q > p
#[(1, 1), (3, 2), (4, 3), (10, 7), (15, 13), (2, 1), (5, 2)]
#clear problemo at r=0.9 ish
#add (7, 3) ind 6
#B^s is large at r=0.6 and r=0.2
#add (8, 5) ind 27 and (16, 15) ind 34
#problemos at the near boundaries
#add (29, 12) ind 37 and (31, 30) ind 44
#significant improvemnet, B^s scale down to 0.002, dJdr scale up to 40 though
#issues now at 0.35 ish and 0.65 ish.
#add (7, 4) ind 19 and (5, 4) ind 41 #both of these should be in earlier.
#now B^s is order 500×10^-6, this is pretty uniform. dJdr is still ~40.
length(curr_surfs)
length(total_surfs)
ininds = [1, 2, 3, 6, 9, 15, 19, 22, 25, 27, 28, 30, 31, 32, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49]
#%%
for i in ininds
    display((total_q[i], total_p[i]))
end
perm = sortperm([total_q[i] / total_p[i] for i in ininds])
rationals = [(total_q[i], total_p[i]) for i in ininds][perm]
#%%
#29 may be ok. (22, 13)
#33 is the same as 39 -> (10, 7)
outinds = [4, 5, 10, 11, 12, 13, 14, 16, 17, 18, 20, 21, 23, 24, 26, 29, 33]
first_surfs = [total_surfs[1], total_surfs[31], total_surfs[40], total_surfs[39], total_surfs[35], total_surfs[3], total_surfs[2], total_surfs[6], total_surfs[27], total_surfs[34], total_surfs[37], total_surfs[44], total_surfs[19], total_surfs[41]];
#first_surfs = [total_surfs[1], total_surfs[31], total_surfs[40], total_surfs[35], total_surfs[3], total_surfs[2], total_surfs[6], total_surfs[27], total_surfs[34], total_surfs[37], total_surfs[44], total_surfs[19], total_surfs[41]];
#now investigate adding more surfaces to non-chaotic region, as it seems possible to get that down to pretty much zero.
#we add [(5, 4), (11, 10), (6, 5), (5, 3), (25, 13), (15, 7), (11, 6), (9, 4)] #again, some of these should probably be in the first pass.
second_surfs = [total_surfs[36], total_surfs[42], total_surfs[43], total_surfs[28], total_surfs[22], total_surfs[15], total_surfs[25], total_surfs[9]];
#with the second surfs, B^s range is ~400×10-6, and is essentially zero other than inside the chaotic region. and dJdr is still ~40.
#now we consider adding more to the chaotic region
#[(11, 8), (17, 11), (23, 14), (13, 10)]
#now these surfaces are actually helping. Starting to see that B^s is v close to zero other than in the islands.
#think this probably does show us that a few more surfaces are needed!
#now add a few extras (7, 5) ind 47 and (16, 13) ind 48, and (14, 11) ind 19
#now it looks like we are at a good spot.
#B^s^2 is basically zero everywhere, except at the islands, which is ~600×10^-9
#dJdr is bigger, but peaks at ~250, which is occuring across the sepratrix at r=0.425 ish and at r=0.575 ish. We may perhaps wish to remove these surfaces.
third_surfs = [total_surfs[38], total_surfs[32], total_surfs[30], total_surfs[45], total_surfs[46], total_surfs[47], total_surfs[48], total_surfs[49]];
#removing each individual surface is still not better
#adding just 30 (23, 14) doesn't make it worse, maybe B^s is a tiny bit better
third_surfs = [total_surfs[38], total_surfs[39]];
#third surfs appear to be determinetal, B^s is pretty much the same, but now dJds is ~200
curr_surfs = total_surfs;
curr_surfs = vcat(surfs4, extra_surfs, surfsc, surfsc2, surfsc3);
curr_surfs = vcat(first_surfs, second_surfs, third_surfs);
curr_surfs = vcat(first_surfs, second_surfs);
#it would be v nice if this could be labeled with the q,p values
#this will probably require changing the surface struct to store p and q. This is a good idea, but perhaps not for now.
#kind of shows that surfs5 is just putting heaps in the 0.7<r<0.9 range, not very useful.
plot_surfs(curr_surfs);
plot_surfs(total_surfs);
B, jac, djac = compute_flux(prob, grids, curr_surfs);
#%%
contourf(θgrid_plot, rgrid_plot, jac[:, :, 1], levels=50, title="Jacobian")
contourf(θgrid_plot, rgrid_plot, djac[1, :, :, 1], levels=50, title="dJdr")
contourf(θgrid_plot, rgrid_plot, djac[2, :, :, 1], levels=50, title="dJdt")
#contourf(θgrid_plot, rgrid_plot, B[1, :, :, 1], levels=50, title="B^s")
contourf(θgrid_plot, rgrid_plot, (B[1, :, :, 1]).^2, levels=50, title="B^s^2")
#%%
Bmean = zeros(Nr);
Bmax = zeros(Nr);
Bsum = zeros(Nr);
djacrmean = zeros(Nr);
djacrmax = zeros(Nr);
djacrsum = zeros(Nr);
djacθmean = zeros(Nr);
djacθmax = zeros(Nr);
djacθsum = zeros(Nr);
for i in 1:Nr
    Bmean[i] = mean((B[1, i, :, :]) .^2)
    Bmax[i] = maximum(B[1, i, :, :])
    Bsum[i] = sum((B[1, i, :, :]) .^2)
    djacrmean[i] = mean((djac[1, i, :, :]) .^2)
    djacrmax[i] = maximum(djac[1, i, :, :])
    djacrsum[i] = sum((djac[1, i, :, :]) .^2)
    djacθmean[i] = mean((djac[2, i, :, :]) .^2)
    djacθmax[i] = maximum(djac[2, i, :, :])
    djacθsum[i] = sum((djac[2, i, :, :]) .^2)
end
#%%
plot(rgrid_plot, Bmean, title="B")
plot(rgrid_plot, Bmax, title="Bmax")
plot(rgrid_plot, Bsum, title="Bsum")
plot(rgrid_plot, djacrmean, title="dJdrmean")
plot(rgrid_plot, djacrmax, title="dJdrmax")
plot(rgrid_plot, djacrsum, title="dJdrsum")
plot(rgrid_plot, djacθmean, title="dJdtmean")
plot(rgrid_plot, djacθmax, title="dJdtmax")
plot(rgrid_plot, djacθsum, title="dJdtsum")
