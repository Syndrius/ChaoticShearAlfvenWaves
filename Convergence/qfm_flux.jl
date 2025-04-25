
#this is essentially testing to see how (B^s)^2 and the qfm jacobian change with the chosen syrfaces, so that we can get an idea of which surfaces to use.
#ideally, we will get some kind of plot showing the values we want, however, paper will probably just outline the process of picking a few, looking at (b^s)^2, adding surfaces where it is large, and repeat, until dJ start to get big.

#need to start being consistent with q/p or p/q etc.
#may be able to replace some of these with smaller values.
#these are in the form a/b, where a is the toroidal winding number, and b is the poloidal winding numebr
#set so rational surface occurs where q(r)=a/b
#17, 11 takes a while to find.
rationals = [(1, 1), (31, 30), (21, 20), (16, 15), (11, 10), (15, 13), (6, 5), (16, 13), (5, 4), (14, 11), (13, 10), (4, 3), (11, 8), (7, 5), (10, 7), (13, 9), (3, 2), (17, 11), (8, 5), (23, 14), (5, 3), (7, 4), (11, 6), (25, 13), (2, 1), (15, 7), (9, 4), (7, 3), (29, 12), (5, 2)]
#ok so I think the surfaces we have are good, perhaps just need to remove some of the surfaces close to the sepratrix, need to check spectrum
#ideally, we want to have a plot of B^2 avg vs number of surfaces, and dJdr vs surfaces and show that there is a sweet spot.
#we just need to figure out how to split up the surfaces.

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
qlocs = [sqrt(i[1]/i[2] - 239/240)/sqrt(5/3) for i in rationals]
#%%
ind = 4
rationals[ind]
qlocs[ind]
for i in 1:length(rationals)
    display(i)
    display(rationals[i])
end
#%%

qlist = [i[1] for i in rationals]
plist = [i[2] for i in rationals]

#%%
total_surfs = construct_surfaces(plist, qlist, qlocs, prob);
save_object("total_surfs.jld2", total_surfs)
#%%
plot_surfs(total_surfs)
#%%
#looks like this is going to require some very specific surface choices.
#start with ~10 surfaces that are small rationals and evenly spread.
first_inds = [1, 5, 7, 12, 14, 17, 21, 25, 28, 30]
for i in first_inds
    display(rationals[i])
end
first_surfs = [total_surfs[i] for i in first_inds];
#%%

#perhaps this would be better, if we just took every 4th for first inds, then every second for this?
#although we want to make the point that not all surfaces are equal.
#this is probably too big of a jump tbh. We may only want to add 3 or 4 surfaces at each iteration.
second_inds = [3, 6, 9, 15, 19, 23, 27]
second_inds = [3, 9, 15, 27]
for i in second_inds
    display(rationals[i])
end
second_surfs = [total_surfs[i] for i in second_inds];
#%%
third_inds = [6, 19, 23]
for i in third_inds
    display(rationals[i])
end
third_surfs = [total_surfs[i] for i in third_inds];
#%%
fourth_inds = [2, 10, 22, 26]
for i in fourth_inds
    display(rationals[i])
end
fourth_surfs = [total_surfs[i] for i in fourth_inds];
#%%
#this is just adding extras away from the chaotic region.
fifth_inds = [4, 8, 24, 29]
for i in fifth_inds
    display(rationals[i])
end
fifth_surfs = [total_surfs[i] for i in fifth_inds];
#%%
#add the more chill ones in the chaotic region
sixth_inds = [13, 16, 20]
for i in sixth_inds
    display(rationals[i])
end
sixth_surfs = [total_surfs[i] for i in sixth_inds];
#%%
#then we can just use the total surfs.
#arr = [1/1, 26/25, 25/23, 11/10, 7/6, 6/5, 5/4, 4/3, 11/8, 7/5, 10/7, 3/2, 8/5, 5/3, 7/4, 9/5, 11/6, 2/1, 13/6, 11/5, 9/4, 7/3, 12/5, 5/2]
#%%
#this function should perhaps be inside QFM.
#probably combine with the compute jac function
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

grids = init_grids(rgrid, θgrid, ζgrid)
rgrid_plot, θgrid_plot, _ = MID.Structures.inst_grids(grids);
#%%
curr_surfs = total_surfs;
curr_surfs = first_surfs;
curr_surfs = vcat(first_surfs, second_surfs);
curr_surfs = vcat(first_surfs, second_surfs, third_surfs);
curr_surfs = vcat(first_surfs, second_surfs, third_surfs, fourth_surfs);
curr_surfs = vcat(first_surfs, second_surfs, third_surfs, fourth_surfs, fifth_surfs);
#total surfs is better than up to sixth lol
#may need to add some more red hot ones to show the Jacobian getting terrible.
#this probably shows the right idea though!
#we will just need to play around with the indices, also need to figure out which surfaces we can actually solve with
#in particular, on Gadi, as mumps seems to be more sensitive.
curr_surfs = vcat(first_surfs, second_surfs, third_surfs, fourth_surfs, fifth_surfs, sixth_surfs);
length(curr_surfs)
curr_surfs = total_surfs;
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
#Bsum = zeros(Nr);
djacrmean = zeros(Nr);
djacrmax = zeros(Nr);
#djacrsum = zeros(Nr);
djacθmean = zeros(Nr);
djacθmax = zeros(Nr);
#djacθsum = zeros(Nr);
for i in 1:Nr
    Bmean[i] = mean((B[1, i, :, :]) .^2)
    Bmax[i] = maximum((B[1, i, :, :]).^2)
    #Bsum[i] = sum((B[1, i, :, :]) .^2)
    djacrmean[i] = mean(abs.(djac[1, i, :, :]))
    djacrmax[i] = maximum(abs.(djac[1, i, :, :]))
    #djacrsum[i] = sum((djac[1, i, :, :]) .^2)
    djacθmean[i] = mean(abs.(djac[2, i, :, :]))
    djacθmax[i] = maximum(abs.(djac[2, i, :, :]))
    #djacθsum[i] = sum((djac[2, i, :, :]) .^2)
end
#%%
plot(rgrid_plot, Bmean, title="B")
plot(rgrid_plot, Bmax, title="Bmax")
#plot(rgrid_plot, Bsum, title="Bsum")
plot(rgrid_plot, djacrmean, title="dJdrmean")
plot(rgrid_plot, djacrmax, title="dJdrmax")
#plot(rgrid_plot, djacrsum, title="dJdrsum")
plot(rgrid_plot, djacθmean, title="dJdtmean")
plot(rgrid_plot, djacθmax, title="dJdtmax")
#plot(rgrid_plot, djacθsum, title="dJdtsum")
