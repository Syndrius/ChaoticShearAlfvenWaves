
#using qfm surfaces for a single island to see what happens
using MID
using MIDViz
using Plots; plotlyjs()
#%%

geo = init_geo(R0=10.0)
k = 0.002
k1 = k
k2 = k

isl = init_island(m0=2, n0=-1, A=k1/2)
#start with no islands
prob = init_problem(geo=geo, q=MID.Equilibrium.island_mode_21, met=:cylinder, isl=isl)
#%%
Ntraj = 150;
rlist = collect(LinRange(0.005, 0.995, Ntraj));
#rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

poincare_plot(prob, Nlaps, Ntraj, rlist)
#%%
rats1 = lowest_rationals(2, prob.q(0.8)[1], prob.q(1.0)[1])
gl1 = surface_guess(rats1, prob.q)
surfs1 = construct_surfaces(rats1, gl1, prob, M=32, N=16);
plot_surfs(surfs1)
#%%
rats2 = lowest_rationals(5, prob.q(0.6)[1], prob.q(0.8)[1])
gl2 = surface_guess(rats2, prob.q)
surfs2 = construct_surfaces(rats2, gl2, prob, M=32, N=16);
plot_surfs(surfs2)
#%%
rats3 = lowest_rationals(7, prob.q(0.4)[1], prob.q(0.6)[1])
gl3 = surface_guess(rats3, prob.q)
surfs3 = construct_surfaces(rats3, gl3, prob, M=32, N=16);
plot_surfs(surfs3)
#%%
rats4 = lowest_rationals(11, prob.q(0.2)[1], prob.q(0.4)[1])
gl4 = surface_guess(rats4, prob.q)
surfs4 = construct_surfaces(rats4, gl4, prob, M=32, N=16);
plot_surfs(surfs4)
#%%
rats5 = lowest_rationals(23, prob.q(0.0)[1], prob.q(0.2)[1])
gl5 = surface_guess(rats5, prob.q)
surfs5 = construct_surfaces(rats5, gl5, prob, M=32, N=16);
plot_surfs(surfs5)
#%%
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5);
#save_object("less_island_surfs.jld2", curr_surfs)
plot_surfs(curr_surfs)
#%%
surfs = load_object("/Users/matt/phd/MID/data/surfaces/island_surfs.jld2");
#this does make things look much better, we may want to be a bit more thorough with this though!
new_surfs = remove_surfs([(15, 7), (13, 6), (11, 5), (13, 7), (11, 6)], surfs);
curr_surfs = surfs;
curr_surfs = new_surfs;
#%%
rgrid_jac = init_grid(type=:rf, N = 100, start=0.03, stop=0.97)
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
minimum(jac)
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
