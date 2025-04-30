
#seems like our original results where significantly better.
#this case will be much closer to our og q profile, where we seemed to get the best results.

#as this should give it the best chance of resolving nicely.

using MID
using MIDViz
using Plots; plotlyjs()
using Plots; gr()
using Statistics
using LaTeXStrings
#%%

function q_prof(r::Float64)
    #chosen so (4, 3)=0.45, (3, 2)=0.55
    a = 239/240
    b = 5/3
    #chosen for (3, 2) isl at 0.5, (2, 1) isl at 0.7
    #a = 47/48
    #b = 25/12
    return a + b*r^2, 2 * b * r
end

#%%
rgrid = init_grid(type=:rf, N = 80, start=0.05, stop=0.95)
θgrid = init_grid(type=:as, start=1, N = 6)
ζgrid = init_grid(type=:as, start=-3, N = 3)

grids = init_grids(rgrid, θgrid, ζgrid)

#%%

geo = init_geo(R0=4.0)
k = 0.001
isl = init_island(m0=4, n0=-3, A = k/15)
isl2 = init_island(m0=3, n0=-2, A = k/3)
isl3 = init_island(m0=2, n0=-1, A=k/2)
prob = init_problem(geo=geo, q=q_prof, isl=isl, isl2=isl2)
#prob = init_problem(geo=geo, q=q_prof, isl=isl3, isl2=isl2)
unprob = init_problem(geo=geo, q=q_prof)

#%%

solver = init_solver(nev=100, targets=[0.23, 0.3, 0.37], prob=unprob)

#%%

evals_norm, ϕ_norm, ϕft_norm = compute_spectrum(prob=unprob, grids=grids, solver=solver);
#%%

continuum_plot(evals_norm, n=-1)
continuum_plot(evals_norm)
ind_norm = find_ind(evals_norm, 0.24897)

potential_plot(ϕft_norm, grids, ind_norm, label_max=0.5)

#%%

Ntraj = 100;
rlist = collect(LinRange(0.4, 0.6, Ntraj));
#rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

#poincare_plot(qfm_benchmark_q, slab_to_plot, Nlaps, Ntraj, 0, 0.0, 0.0, 0.0, R0, isl, MID.Structures.no_isl, rlist)
poincare_plot(prob, Nlaps, Ntraj, rlist,  prob.geo.R0)#, filename=save_dir * "original_poincare.png");

#%%
q_prof(1.0)
rat1 = [(1, 1), (2, 1), (3, 2), (5, 2), (4, 3), (5, 3), (7, 3), (5, 4), (7, 4), (9, 4), (6, 5), (7, 5), (8, 5), (9, 5), (11, 5), (12, 5), (13, 5)]
ql1 = [i[1] for i in rat1]
pl1 = [i[2] for i in rat1]
g1 = [sqrt(i[1]/i[2] - 239/240) / sqrt(5/3) for i in rat1]
@time surfs1 = construct_surfaces(pl1, ql1, g1, prob);
plot_surfs(surfs1)
#%%
rat2 = [(31, 30), (21, 20), (11, 10)]
ql2 = [i[1] for i in rat2]
pl2 = [i[2] for i in rat2]
g2 = [sqrt(i[1]/i[2] - 239/240) / sqrt(5/3) for i in rat2]
@time surfs2 = construct_surfaces(pl2, ql2, g2, prob);

plot_surfs(surfs2)
#%%
rat3 = [(10, 7), (11, 7), (11, 8), (9, 7)]
ql3 = [i[1] for i in rat3]
pl3 = [i[2] for i in rat3]
g3 = [sqrt(i[1]/i[2] - 239/240) / sqrt(5/3) for i in rat3]
@time surfs3 = construct_surfaces(pl3, ql3, g3, prob);

plot_surfs(surfs3)
#%%
rat4 = [(17, 15), (41, 40)]
ql4 = [i[1] for i in rat4]
pl4 = [i[2] for i in rat4]
g4 = [sqrt(i[1]/i[2] - 239/240) / sqrt(5/3) for i in rat4]
@time surfs4 = construct_surfaces(pl4, ql4, g4, prob);

plot_surfs(surfs4)
#%%
rat5 = [(23, 20)]
ql5 = [i[1] for i in rat5]
pl5 = [i[2] for i in rat5]
g5 = [sqrt(i[1]/i[2] - 239/240) / sqrt(5/3) for i in rat5]
@time surfs5 = construct_surfaces(pl5, ql5, g5, prob);

plot_surfs(surfs5)
#%%
tmp_surfs = vcat(surfs3[1:2], [surfs3[4]]);
tmp_surfs = surfs3[1:2];
curr_surfs = vcat(surfs1, surfs2, tmp_surfs, surfs4, surfs5);
plot_surfs(curr_surfs)
#%%
evals, ϕ, ϕft = compute_spectrum_qfm(grids=grids, prob=prob, solver=solver, surfs=curr_surfs);

continuum_plot(evals, n=-1)
ind = find_ind(evals, 0.294)
potential_plot(ϕft, grids, ind)
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
    #ξr, wgr = MID.Construct.FastGaussQuadrature.gausslegendre(grids.r.gp) #same as python!
    #ξθ, wgθ = MID.Construct.FastGaussQuadrature.gausslegendre(grids.θ.gp)
    #ξζ, wgζ = MID.Construct.FastGaussQuadrature.gausslegendre(grids.ζ.gp)

    #struct for storing the intermediate data for the coordinate transform
    CT = MID.QFM.CoordTsfmT()
    #jac = zeros(length(rvals), length(θvals), length(ζvals))
    jac = zeros(length(rgrid), length(θgrid), length(ζgrid))
    djac = zeros(3, length(rgrid), length(θgrid), length(ζgrid))
    B = zeros(3, length(rgrid), length(θgrid), length(ζgrid))

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
B, jac, djac = compute_flux(prob, grids, curr_surfs);
#%%
rgrid_plot, θgrid_plot, _ = MID.Structures.inst_grids(grids);
contourf(θgrid_plot, rgrid_plot, jac[:, :, 1], levels=50, title="Jacobian")
contourf(θgrid_plot, rgrid_plot, djac[1, :, :, 1], levels=50, title="dJdr")
contourf(θgrid_plot, rgrid_plot, djac[2, :, :, 1], levels=50, title="dJdt")
#contourf(θgrid_plot, rgrid_plot, B[1, :, :, 1], levels=50, title="B^s")
contourf(θgrid_plot, rgrid_plot, (B[1, :, :, 1]).^2, levels=50, title="B^s^2")
#%%
B2mean = zeros(length(surfs_list));
djacdrmean = zeros(length(surfs_list));
djacdθmean = zeros(length(surfs_list));

for i in 1:length(surfs_list)
    curr_surfs = surfs_list[i]

    B, jac, djac = compute_flux(prob, grids, curr_surfs);

    B2mean[i] = mean(B[1, :, :, :] .^ 2)
    djacdrmean[i] = mean(abs.(djac[1, :, :, :]))
    djacdθmean[i] = mean(abs.(djac[2, :, :, :]))
end
#%%
