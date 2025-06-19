
#using qfm surfaces for a single island to see what happens
#think flux is better mainly just for behavioour at the axis.
#and the q-profiel is less cooked
#for the radial case we will probably need to change the amplitude of the island
#which will be annoying!
using MID
using MIDViz
using Statistics
using JLD2
using Plots; gr()
using Plots; plotlyjs()
#%%

geo = init_geo(R0=4.0)

isl21a = init_island(m0=2, n0=-1, w=0.1, ψ0=0.5, qp=2.0, flux=true)
#start with no islands
prob = init_problem(geo=geo, q=island_q, met=:cylinder, isl=isl21a, type=:flux)
#%%
rvals = LinRange(0.0, 1.0, 100)
qvals = [prob.q(i)[1] for i in rvals]
plot(rvals, qvals)
#%%
Ntraj = 150;
rlist = collect(LinRange(0.005, 0.995, Ntraj));
#rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

poincare_plot(prob, Nlaps, Ntraj, rlist)
#%%
rats1 = lowest_rationals(5, prob.q(0.8)[1], prob.q(1.0)[1])
gl1 = surface_guess(rats1, prob.q)
surfs1 = construct_surfaces(rats1, gl1, prob, M=32, N=16);
plot_surfs(surfs1)
#%%
rats2 = lowest_rationals(6, prob.q(0.6)[1], prob.q(0.8)[1])
gl2 = surface_guess(rats2, prob.q)
surfs2 = construct_surfaces(rats2, gl2, prob, M=32, N=16);
plot_surfs(surfs2)
#%%
#this will be the spiciest one!
#we will probably remove all of these
#also need to determine if we want to keep the central one or not.
#I think we do!
#we may want to look at the coordinate maps to see if we can make sure it is mapping as expected.
rats3 = lowest_rationals(7, prob.q(0.4)[1], prob.q(0.6)[1])
gl3 = surface_guess(rats3, prob.q)
surfs3 = construct_surfaces(rats3, gl3, prob, M=32, N=16);
plot_surfs(surfs3)
#%%
rats4 = lowest_rationals(9, prob.q(0.2)[1], prob.q(0.4)[1])
gl4 = surface_guess(rats4, prob.q)
surfs4 = construct_surfaces(rats4, gl4, prob, M=32, N=16);
plot_surfs(surfs4)
#%%
rats5 = lowest_rationals(13, prob.q(0.0)[1], prob.q(0.2)[1])
gl5 = surface_guess(rats5, prob.q)
surfs5 = construct_surfaces(rats5, gl5, prob, M=32, N=16);
plot_surfs(surfs5)
#%%
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5);
save_object("/Users/matt/phd/MID/data/surfaces/flux_island_21b.jld2", curr_surfs)
plot_surfs(curr_surfs)
#%%
surfs = load_object("/Users/matt/phd/MID/data/surfaces/island_surfs.jld2");
#this does make things look much better, we may want to be a bit more thorough with this though!
new_surfs = remove_surfs([(15, 7), (13, 6), (11, 5), (13, 7), (11, 6)], surfs);
curr_surfs = surfs;
curr_surfs = new_surfs;
#%%
#Jacobian is propto R0, so it is a bit hard to tell here if the jacobian is actually to large?
#we may need to try with lower R0, but then matching with Zhisongs code could be q bit annoying?
#I guess we can always just match it with ours?
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
#%%
sgrid = init_grid(type=:rf, N=100, sep1=0.44, sep2=0.56, frac=0.5)
ϑgrid = init_grid(type=:as, N=4, start=1)
φgrid = init_grid(type=:as, N=2, start=-2)
grids = init_grids(sgrid, ϑgrid, φgrid)
#%%
solver = init_solver(prob=prob, full_spectrum=true)
#%%
evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=curr_surfs)

#%%
#hard to tell if spire is due to island or cooked jacobian!
continuum_plot(evals)

island_om = evals.ω[0.44 .< evals.x1 .< 0.56]
#%%
#isl_ind=1 gives a pretty decent island mode.
#so do most of them, looks like the qfm case will at least make the results more tractable.
#given this is a tiny grid.
#also no zero zero harmonic in this test.
#think we have at least the maximim surfaces. May need to remove some around the sepratrix as per.
isl_ind = 145

ind = find_ind(evals, island_om[isl_ind])
ind = find_ind(evals, 0.190136)
potential_plot(ϕft, grids, ind)

contour_plot(ϕ, grids, ind=ind)


#%%
