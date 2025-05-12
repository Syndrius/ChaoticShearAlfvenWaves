#this will hopefully be a more organised way of picking the surfaces
#note the we probably need to fix qfm density!
#ok so trying to use a q-profile that would work for tae and islands with the same n
#is not going to work, even with a density profile
#so, either we create a chaotic region with a variety of islands, 
#or we use a rapidly increasing q-profile, and accepts that the tae will look shit.

using MID
using MIDViz
using JLD2
using Plots; plotlyjs()
using Plots; gr()
#using Roots
using Statistics
#%%
#=
function low_shear_qfm_q(r::Float64)
    a = 1.2
    b = 0.6
    c = 0.1
    return a + b *r^2 + c*r^4, 2*b*r + 4*c*r^3
end
=#
#first we looks at the poincare plot

k1 = 0.02
geo = init_geo(R0=4.0)
isl1 = init_island(m0=3, n0=2, A=k1/3)

isls = [isl1]#, isl2, isl3]

prob = init_problem(geo=geo, q=low_shear_qfm_q, isls=isls)

#%%

Ntraj = 150
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

poincare_plot(prob, Nlaps, Ntraj, rlist)
#%%
#so (11, 7) does not work!, regardless of res
rats1 = lowest_rationals(7, prob.q(0.0)[1], prob.q(1.0)[1])
gl1 = surface_guess(rats1, prob.q)
#changing these numbers doesn't really help remove the spikes
@time surfs1 = construct_surfaces(rats1, gl1, prob, M=32, N=16);
plot_surfs(surfs1)
#%%
#now we need more surfs below 0.25
rats2 = lowest_rationals(25, prob.q(0.0)[1], prob.q(0.25)[1])
gl2 = surface_guess(rats2, prob.q)
surfs2 = construct_surfaces(rats2, gl2, prob, M=32, N=16);
plot_surfs(surfs2)
#%%
#still need more below 0.2
#gross
rats3 = [(53, 44), (13, 9), (14, 9), (13, 8), (13, 10), (14, 11)]
gl3 = surface_guess(rats3, prob.q)
surfs3 = construct_surfaces(rats3, gl3, prob, M=32, N=16);
plot_surfs(surfs3)
#%%
#shouldn't add (14, 11) here
#rats4 = [(11, 8), (14, 11), (26, 17), (22, 15)]
rats4 = [(11, 8), (26, 17), (22, 15)]
gl4 = surface_guess(rats4, prob.q)
surfs4 = construct_surfaces(rats4, gl4, prob, M=32, N=16);
plot_surfs(surfs4)
#%%
rats5 = [(19, 15), (26, 21), (15, 11), (17, 13)]
gl5 = surface_guess(rats5, prob.q)
surfs5 = construct_surfaces(rats5, gl5, prob, M=32, N=16);
plot_surfs(surfs5)
#%%
#add more to chaotic region, most likely these will make it worse.
#these surfs make B^s ~10x larger in chaotic region.
#somewhat expected
rats6 = [(18, 11), (17, 10), (16, 9), (20, 11)]
gl6 = surface_guess(rats6, prob.q)
surfs6 = construct_surfaces(rats6, gl6, prob, M=32, N=16);
plot_surfs(surfs6)
#%%
rats7 = lowest_rationals(60, prob.q(0.0)[1], prob.q(0.1)[1])
#these are a bit stupid but might fix the axis problemos.
#otherwise we can just got from 0.05
rats7 = [(41, 34), (71, 59)]
gl7 = surface_guess(rats7, prob.q)
#changing these numbers doesn't really help remove the spikes
@time surfs7 = construct_surfaces(rats7, gl7, prob, M=32, N=16);
plot_surfs(surfs7)
#%%
curr_surfs = vcat(surfs1, surfs2);
curr_surfs = vcat(surfs1, surfs2, surfs3);
#19/7 surface is overlapping with neighbour.
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4);
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5);
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4[1:1], surfs4[3:end], surfs5);
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5, surfs6);

curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5, surfs6, surfs7);

#(11, 7) is cooke d for some reason. Bit surprising as it looks pretty fine.
#just had it twice lol.
#curr_surfs = vcat(surfs1, surfs2, surfs3[1:1], surfs3[3:6]);
#curr_surfs = vcat(surfs1, surfs2, surfs3[1:1], surfs3[3:6], surfs4);
#curr_surfs = vcat(surfs1, surfs2, surfs3[1:1], surfs3[3:6], surfs4, surfs5);
#curr_surfs = vcat(surfs1, surfs2, surfs3[1:1], surfs3[3:6], surfs4, surfs5, surfs6);
plot_surfs(curr_surfs)

save_object("bench_low_shear_surfs.jld2", curr_surfs)

#%%
rgrid_jac = init_grid(type=:rf, N = 100)#, start=0.0, stop=0.95)
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
rgrid = init_grid(type=:rf, N = 150)#, sep1=0.7, sep2=0.9, frac=0.4)
θgrid = init_grid(type=:as, N = 11, start=1) 
ζgrid = init_grid(type=:as, N = 3, start=-3)
grids = init_grids(rgrid, θgrid, ζgrid)
solver = init_solver(nev=150, targets=[0.20, 0.25, 0.30, 0.35], prob=prob)
evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=curr_surfs);

#%%

continuum_plot(evals, legend=false)#, n=-2)
ind = find_ind(evals, 0.2866)
potential_plot(ϕft, grids, ind, label_max=0.05)
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
