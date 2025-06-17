#we will do this in flux coordinates as well.

#testing the cylindrical stuff in radial coordinates
#im sure it will def work now 
#yep these surfs are also completly cooked
#very nice
#Maybe something is wrong with qfm?
#don't think there are any actual examples of qfm being used in completly chaotic regions
#we may have to try with small chaotic regions, and just ignore the peices where the islands still exist.
#which will be quite shite.
#so perhaps our previous example may have been ok.
#think we will need to just consider the more regular case of a wee bit chaotic.
#very shit tbh.


using MID
using MIDViz
using JLD2
using Plots; plotlyjs()
using Plots; gr()
using Statistics
#%%
#first we looks at the poincare plot
function cyl_qfm_q(r::Float64)
    return 1.0 + r, 1.0
end
#%%
k1 = 0.00025
#k2 = 0.0003
#k3 = 0.0000
geo = init_geo(R0=4.0)
isl1 = init_island(m0=3, n0=-2, A=k1)
isl2 = init_island(m0=5, n0=-3, A=k1)
#isl3 = init_island(m0=8, n0=-5, A=k3/8) #unsure if we will want this one as well

isls = MID.Geometry.IslandT[isl1, isl2]

prob = init_problem(geo=geo, q=cyl_qfm_q, isls=isls, met=:cylinder)

#%%

Ntraj = 300;
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
rlist = collect(LinRange(0.4, 0.8, Ntraj));
Nlaps = 500;

poincare_plot(prob, Nlaps, Ntraj, rlist, ylimits=(0.5, 0.7))
#%%
#so (11, 7) does not work!, regardless of res
rats1 = lowest_rationals(11, prob.q(0.0)[1], prob.q(1.0)[1])
gl1 = surface_guess(rats1, prob.q)
#changing these numbers doesn't really help remove the spikes
surfs1 = construct_surfaces(rats1, gl1, prob, M=32, N=16);
plot_surfs(surfs1)
#%%
#now we need more surfs below 0.05 and above 0.9 ish
#seems to be taking years for some reason
#ok so maybe we just dont worry about these.
#so 20, 19 is just a bit cooked for some reason!
#30 29 is also not findable..
#surfaces seem to get a bit mroe cooked at r=0...
#probs because we don't have any island amp modification.
rats2 = [(16, 15), (29, 15), (39, 20)]
#rats2 = [(30, 29), (23, 12), (39, 20)]
#rats2 = [(23, 12), (16, 15), (29, 15)]
#rats2 = [(39, 20)]
gl2 = surface_guess(rats2, prob.q)
surfs2 = construct_surfaces(rats2, gl2, prob, M=32, N=16);
plot_surfs(surfs2)
#%%
rats2 = lowest_rationals(16, prob.q(0.91)[1], prob.q(1.0)[1])
curr_surfs = surfs1;
curr_surfs = vcat(surfs1, surfs2);
curr_surfs = vcat(surfs1, surfs2, surfs3);
#19/7 surface is overlapping with neighbour.
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4);
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5);
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4[1:1], surfs4[3:end], surfs5);
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5, surfs6, surfs7);
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5, surfs7);

#(11, 7) is cooke d for some reason. Bit surprising as it looks pretty fine.
#just had it twice lol.
#curr_surfs = vcat(surfs1, surfs2, surfs3[1:1], surfs3[3:6]);
#curr_surfs = vcat(surfs1, surfs2, surfs3[1:1], surfs3[3:6], surfs4);
#curr_surfs = vcat(surfs1, surfs2, surfs3[1:1], surfs3[3:6], surfs4, surfs5);
#curr_surfs = vcat(surfs1, surfs2, surfs3[1:1], surfs3[3:6], surfs4, surfs5, surfs6);
plot_surfs(curr_surfs)

#save_object("low_shear_surfs.jld2", curr_surfs)
curr_surfs = load_object("low_shear_surfs.jld2");
save_object("/Users/matt/phd/MID/data/cyl_qfm_surfaces/rad_surfaces.jld2", curr_surfs)

#%%
rgrid_jac = init_grid(type=:rf, N = 100, start=0.1, stop=0.9)
θgrid_jac = init_grid(type=:af, N = 20) 
ζgrid_jac = init_grid(type=:af, N = 4)
grids_jac = init_grids(rgrid_jac, θgrid_jac, ζgrid_jac)
B, jac, djac = compute_jac(prob, grids_jac, curr_surfs);
#%%
rgrid_plot, θgrid_plot, ζgrid_plot = MID.Structures.inst_grids(grids_jac);

contourf(θgrid_plot, rgrid_plot, jac[:, :, 1], levels=50, title="Jacobian")
contourf(θgrid_plot, rgrid_plot, djac[1, :, :, 1], levels=50, title="dJdr")
contourf(θgrid_plot, rgrid_plot, djac[2, :, :, 1], levels=50, title="dJdt")
#contourf(θgrid_plot, rgrid_plot, B[1, :, :, 1], levels=50, title="B^s")
contourf(θgrid_plot, rgrid_plot, (B[1, :, :, 1]).^2, levels=50, title="B^s^2")
#%%
Nr = rgrid_jac.N
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

#%%
rgrid = init_grid(type=:rf, N = 150, start=0.1, stop=0.9)
θgrid = init_grid(type=:as, N = 11, start=1) 
ζgrid = init_grid(type=:as, N = 3, start=-3)
grids = init_grids(rgrid, θgrid, ζgrid)
solver = init_solver(nev=150, targets=[0.20, 0.25, 0.30, 0.35], prob=prob)
evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=curr_surfs);

#%%

#this probably already shows adequate results.
#we will just have to ignore the islands etc.
#but if we zoom in very closely to the chaotic region
#we get normal continuum modes.
continuum_plot(evals, legend=false)#, n=-2)
ind = find_ind(evals, 0.2756605)
potential_plot(ϕft, grids, ind, label_max=0.05)
contour_plot(ϕ, grids, ind=ind, label_max=0.05)
#%%

unprob = init_problem(geo=geo, q=cyl_qfm_q, met=:cylinder)
unevals, unϕ, unϕft = compute_spectrum(prob=unprob, grids=grids, solver=solver);
continuum_plot(unevals, legend=false)#, n=-2)

#%%
using Plots; gr()

Ntraj = 200;
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
rlist = collect(LinRange(0.4, 0.8, Ntraj));
Nlaps = 500;

x, z = poincare_plot(prob, Nlaps, Ntraj, rlist, ylimits=(0.5, 0.7))

MIDViz.plot_surface_poincare(x, z, curr_surfs, ylimits=(0.5, 0.66))

#%%
plot_font = "Computer Modern"
default(fontfamily=plot_font, grid=false, framestyle=:semi, palette=:tol_bright)
cur_colors = palette(:tol_bright);
lfs = 18
xfs = 22
yfs = 22
tfs = 16
ylim= (0.48, 0.7)
MIDViz.plot_poincare(x, z, color=:black, xguidefontsize=xfs, yguidefontsize=yfs, xtickfontsize=tfs, ytickfontsize=tfs-2, dpi=1200, ylimits=ylim)
MIDViz.Plotting.overlay_surfs(curr_surfs, color=cur_colors[2], linewidth=2.5)

savefig("aapps_abstract_poincare.png")
