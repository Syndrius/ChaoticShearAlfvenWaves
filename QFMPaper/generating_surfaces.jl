
#so this just straight up doesn't work
#nice, real nice
using MID
using MIDViz
using JLD2
using Plots; gr()
using Plots; plotlyjs()
#%%
function cyl_qfm_q(ψ::Float64)

    return 1.0 + ψ, 1.0
end
#%%
geo = init_geo(R0=10.0)
isl1 = init_island(m0=7, n0=-5, A=0.006/7)
isl2 = init_island(m0=3, n0=-2, A=0.004/3)
isl3 = init_island(m0=8, n0=-5, A=0.006/8)

#isl1 = init_island(m0=7, n0=-5, A=0.001/7)
#isl2 = init_island(m0=3, n0=-2, A=0.001/3)
#isl3 = init_island(m0=8, n0=-5, A=0.001/8)
isls = [isl1, isl2, isl3]

prob = init_problem(geo=geo, q=cyl_qfm_q, isls=isls, met=:cylinder)
#%%
Ntraj = 100
rlist = collect(LinRange(0.05, 0.95, Ntraj));
#rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 800;

x, z = poincare_plot(prob, Nlaps, Ntraj, rlist)#, color=:black)
#%%
rats1 = lowest_rationals(13, prob.q(0.0)[1], prob.q(1.0)[1])
#rats1 = [(8, 5)]
#rats1 = vcat(rats1[1:6], rats1[8:end])
gl1 = surface_guess(rats1, prob.q)
#gl1 = 0.5 .* ones(length(rats1))
#changing these numbers doesn't really help remove the spikes
surfs1 = construct_surfaces(rats1, gl1, prob, M=32, N=32);
plot_surfs(surfs1)
#%%
#need more at edges
rats2 = [(21, 20), (31, 30), (29, 15), (39, 20)]
#rats1 = [(8, 5)]
#rats1 = vcat(rats1[1:6], rats1[8:end])
gl2 = surface_guess(rats2, prob.q)
#gl1 = 0.5 .* ones(length(rats1))
#changing these numbers doesn't really help remove the spikes
surfs2 = construct_surfaces(rats2, gl2, prob, M=32, N=16);
plot_surfs(surfs2)
#%%
#need more at edges
rats3 = [(14, 9), (13, 9)]
#rats1 = [(8, 5)]
#rats1 = vcat(rats1[1:6], rats1[8:end])
gl3 = surface_guess(rats3, prob.q)
#gl1 = 0.5 .* ones(length(rats1))
#changing these numbers doesn't really help remove the spikes
surfs3 = construct_surfaces(rats3, gl3, prob, M=64, N=8);
plot_surfs(surfs3)
#%%
curr_surfs = vcat(surfs1, surfs2);
plot_surfs(curr_surfs)
#%%
#probably as good as we can hope for.
to_remove = [(18, 13), (21, 13), (17, 12), (19, 12), (20, 13), (19, 13)]
new_surfs = MID.QFM.QFMSurfaceT[]
for surf in curr_surfs
    if !(surf.rational in to_remove)
        push!(new_surfs, surf)
    end
end
curr_surfs = new_surfs;
plot_surfs(curr_surfs)
#%%
rgrid_jac = init_grid(type=:rf, N = 100, start=0.05, stop=0.95)
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
rgrid = init_grid(type=:rf, N = 150, start=0.05, stop=0.95)#, sep1=0.7, sep2=0.9, frac=0.4)
θgrid = init_grid(type=:as, N = 11, start=1) 
ζgrid = init_grid(type=:as, N = 3, start=-3)
grids = init_grids(rgrid, θgrid, ζgrid)
solver = init_solver(nev=150, targets=[0.20, 0.25, 0.30, 0.35], prob=prob)
evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=curr_surfs);

#%%

continuum_plot(evals, legend=false)#, n=-2)
ind = find_ind(evals, 0.32240)
potential_plot(ϕft, grids, ind, label_max=0.05)
contour_plot(ϕ, grids, ind=ind)
#%%

unprob = init_problem(q=prob.q, dens=dens_prof, geo=geo)

evals_norm, ϕ_norm, ϕft_norm = compute_spectrum(prob=unprob, grids=grids, solver=solver);
#%%
continuum_plot(evals_norm, legend=false)#, n=-2)

ind_norm = find_ind(evals_norm, 0.2327)

potential_plot(ϕft_norm, grids, ind_norm, label_max=0.2)

#%%


rgrid_cont = init_grid(type=:rc, N=150)
θgrid_cont = init_grid(type=:as, N = 7, start=1) 
ζgrid_cont = init_grid(type=:as, N = 3, start=-3)
grids_cont = init_grids(rgrid_cont, θgrid_cont, ζgrid_cont)

evals_cont = compute_continuum(unprob, grids_cont);

#%%

continuum_plot(evals_cont, grids_cont)
#%%

B2mean = zeros(length(surfs_list));
djacdrmean = zeros(length(surfs_list));
djacdθmean = zeros(length(surfs_list));

for i in 1:length(surfs_list)
    display(i)
    curr_surfs = surfs_list[i]

    B, jac, djac = compute_jac(prob, grids_jac, curr_surfs);

    B2mean[i] = mean(B[1, :, :, :] .^ 2)
    djacdrmean[i] = mean(abs.(djac[1, :, :, :]))
    djacdθmean[i] = mean(abs.(djac[2, :, :, :]))
end
