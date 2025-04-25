#here we offer some tests for properly choosing the right number and specific surfaces.
#this is not very scientific, ideally we would have more.

#this will mainly consist of looking at the new poincare plot
#and looking at the jacobain in the new coordinates, where we want the value to not get that large.

#ideally, in the future, we will have a much more sophisticated method of choosing the surfaces.
using MID
using MIDViz
using Plots
using JLD2
#%%
R0=4.0

geo = init_geo(R0=R0)

#non-resonant case
k = 0.05
isl = init_island(m0=3, n0=2, A=k/5)
isl2 = init_island(m0=7, n0=-3, A=0.0)
prob = init_problem(q=qfm_q, geo=geo, isl=isl, isl2=isl2)#, flr=flr)

#chaotic case
k = 0.0006
isl = init_island(m0=3, n0=-2, A=k/3)
isl2 = init_island(m0=4, n0=-3, A=k/4)
#isl2 = init_island(m0=4, n0=-3, A=0.0)

prob = init_problem(q=qfm_q, geo=geo, isl=isl, isl2=isl2)#, flr=flr)
#%%
qlist, plist = farey_tree(5, 1, 1, 2, 1)
guess_list = 0.5 .* ones(length(qlist));
guess_list = @. sqrt(qlist / plist - 0.95) / sqrt(1.5);
#needed for the single island chain case, this one is probably at the sepratrix I guess.
@time surfs = construct_surfaces(plist, qlist, guess_list, prob);
plot_surfs(surfs)
#%%
#save_object("surf_test_chaos.jld2", surfs)
save_object("surf_test_isl.jld2", surfs)
#these are worse based on Jacobian eye test. Ideally we could make this a bit more analytical.
save_object("surf_test_chaos5.jld2", surfs)
#%%
surfs_chaos = load_object("surf_test_chaos.jld2");
surfs_bench = load_object("surf_test_bench.jld2");
plot_surfs(surfs_chaos)
#%%
#the 20/21 surface has a huge effect on the jacobian around 0.2, before it was varying by like 2
#now it varies by 0.2.
#so this doesn't work for q/p < 1 or whichever way around.
#poses a problem for the current q-profile, where 1/1 is at 0.2 ish, meaning we cannot put another surface at 0.05 etc, which would be required to get the full domain.
#the 12/29 surfaces is at 0.97 ish, this fixes the problems for r>0.8.
#looks like we will need a new q-profile where the 1,1 surface is at 0.05 ish.
extra_surfs = construct_surfaces([15, 13, 20, 12], [16, 15, 21, 29], [0.3, 0.3, 0.26, 0.95], prob);
comb_surfs = vcat(extra_surfs, surfs_chaos);
plot_surfs(comb_surfs)
#%%
#think we should actually put this inside MID, as this is quite a useful function for QFM checking!
function compute_jac(prob, grids, surfs)

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
    jac_tor = zeros(grids.r.N, grids.θ.N, grids.ζ.N)
    djac = zeros(3, grids.r.N, grids.θ.N, grids.ζ.N)

    #for (i, r) in enumerate(rvals), (j, θ) in enumerate(θvals), (k, ζ) in enumerate(ζvals)
    for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)
        MID.QFM.coord_transform!(r, θ, ζ, CT, surf_itp, sd)
        MID.Geometry.toroidal_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], 4.0)
        MID.Equilibrium.compute_B!(tor_B, tor_met, prob.q, prob.isl, prob.isl2, CT.coords[1], CT.coords[2], CT.coords[3])
        MID.QFM.met_transform!(tor_met, qfm_met, CT)
        MID.QFM.B_transform!(tor_B, qfm_B, qfm_met, CT)

        jac[i, j, k] = qfm_met.J[1]
        djac[:, i, j, k] = qfm_met.dJ[:]
        jac_tor[i, j, k] = tor_met.J[1]
        #jac[i, j, k] = qfm_B.B[1]
        #jac_tor[i, j, k] = qfm_B.B[2]
    end
    return jac, djac, jac_tor
end
#%%



Nr = 100
Nθ = 20
Nζ = 1
#ok so the jacobian below ~r=0.4 is completly cooked. No idea why. unsure how to fix.
#I guess this is because of the axis? not really sure how to fix that though?
rgrid = init_grid(type=:rf, N = Nr, start = 0.2, stop =0.9999)
θgrid = init_grid(type=:af, N = Nθ, pf=5)
ζgrid = init_grid(type=:af, N = Nζ, pf=-2)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%
jac, djac, jac_tor = compute_jac(prob, grids, comb_surfs);
jac, jac_tor = compute_jac(prob, grids, surfs);
jac, jac_tor = compute_jac(prob, grids, surfs_chaos);
#%%
rgrid_plot, θgrid_plot, _ = MID.Structures.inst_grids(grids);
contourf(θgrid_plot, rgrid_plot, jac[:, :, 1], levels=50)
contourf(θgrid_plot, rgrid_plot, djac[1, :, :, 1], levels=50)
contourf(θgrid_plot, rgrid_plot, djac[2, :, :, 1], levels=50)
contourf(θgrid_plot, rgrid_plot, jac_tor[:, :, 1], levels=50)
#%%
Ntraj = 100
rlist = collect(LinRange(0.0001, 1.0, Ntraj));
Nlaps = 500
poincare_plot(prob, Nlaps, Ntraj, rlist, R0);
#%%
#now the qfm poincare plot
Ntraj = 40
#starting this below 0.4 causes immediate issues for the Jacobain being enourmous/negative
rlist = collect(LinRange(0.4, 0.9, Ntraj));
Nlaps = 500
x, z = poincare_plot(prob, Nlaps, Ntraj, rlist, R0, surfs=surfs_chaos);
#%%
#now the qfm poincare plot with better surfaces
Ntraj = 40
#this doesn't seem to be any better.
#I guess this doesn't matter if the range has to be restricted tho!
rlist = collect(LinRange(0.4, 0.9, Ntraj));
Nlaps = 500
x, z = poincare_plot(prob, Nlaps, Ntraj, rlist, R0, surfs=comb_surfs);



#%%

qlist = zeros(100);
rlist = LinRange(0, 0.5, 100);
for (i, r) in enumerate(rlist)
    qlist[i], _ = qfm_q(r)
end
#%%

plot(rlist, qlist)
