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
using Roots
using Statistics
#%%
function q_prof(r::Float64)
    a = 1.2
    b = 0.6
    c = 0.1
    d = 0.0
    return a + b*r^2 + c*r^4 + d*r^6, 2 * b * r + 4*c*r^3 + 6*d*r^5
    #return a + b*r^2 + c*r^10 + d*r^4, 2*b*r + 10*c*r^9 + 4*d*r^3
    #return a + b*r^6, 6 * b * r^5
    #return a + b*r^4, 4 * b * r^3
end

function dens_prof(r::Float64)
    a = 1.05
    b = 0.4
    c = 5.0
    d = 6

    #0.8 and 0.2 give a reasnable approximate to the density decrease at the edge
    #while also maintaing the tae gap throughout the domain
    #note that this does not prevent continuum interactions.
    rh = 0.8
    scale = 0.2

    return 1.0 #1/2 * (1-tanh((r-rh)/scale))

end
function surf_guess(rats, q_prof)
    gl = Float64[]

    for i in rats
        diff(r) = i[1]/i[2] - q_prof(r)[1]
        sol = find_zero(diff, 0.5)
        push!(gl, sol)
    end
    return gl
end


lowest_rationals(20, q_prof(0.64)[1], q_prof(0.72)[1])
#%%
#first we looks at the poincare plot

k1 = 0.0019
k2 = 0.0007
k3 = 0.0013
geo = init_geo(R0=4.0)
isl1 = init_island(m0=7, n0=-4, A=k1/7)
isl2 = init_island(m0=5, n0=-3, A=k2/5)
isl3 = init_island(m0=8, n0=-5, A=k3/8) #unsure if we will want this one as well

isls = [isl1, isl2, isl3]

prob = init_problem(geo=geo, q=q_prof, dens=dens_prof, isls=isls)

#%%

Ntraj = 150;
rlist = collect(LinRange(0.7, 0.95, Ntraj));
#rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

poincare_plot(prob, Nlaps, Ntraj, rlist)
#%%
#so (11, 7) does not work!, regardless of res
rats1 = lowest_rationals(7, q_prof(0.0)[1], q_prof(1.0)[1])
alist1 = [i[1] for i in rats1]
blist1 = [i[2] for i in rats1]
gl1 = surf_guess(rats1, q_prof)
#changing these numbers doesn't really help remove the spikes
surfs1 = MID.QFM.new_construct_surfaces(alist1, blist1, gl1, prob, M=32, N=16);
plot_surfs(surfs1)
#%%
#now we need more surfs below 0.4
rats2 = lowest_rationals(25, q_prof(0.0)[1], q_prof(0.25)[1])
alist2 = [i[1] for i in rats2]
blist2 = [i[2] for i in rats2]
gl2 = surf_guess(rats2, q_prof)
surfs2 = MID.QFM.new_construct_surfaces(alist2, blist2, gl2, prob, M=32, N=16);
plot_surfs(surfs2)
#%%
#still need more below 0.2
#gross
rats3 = [(53, 44), (11, 7), (13, 9), (14, 9), (13, 8), (13, 10), (14, 11)]
alist3 = [i[1] for i in rats3]
blist3 = [i[2] for i in rats3]
gl3 = surf_guess(rats3, q_prof)
surfs3 = MID.QFM.new_construct_surfaces(alist3, blist3, gl3, prob, M=32, N=16);
plot_surfs(surfs3)
#%%
rats4 = [(11, 8), (14, 11), (26, 17), (22, 15)]
alist4 = [i[1] for i in rats4]
blist4 = [i[2] for i in rats4]
gl4 = surf_guess(rats4, q_prof)
surfs4 = MID.QFM.new_construct_surfaces(alist4, blist4, gl4, prob, M=32, N=16);
plot_surfs(surfs4)
#%%
rats5 = [(19, 15), (26, 21), (15, 11), (17, 13)]
alist5 = [i[1] for i in rats5]
blist5 = [i[2] for i in rats5]
gl5 = surf_guess(rats5, q_prof)
surfs5 = MID.QFM.new_construct_surfaces(alist5, blist5, gl5, prob, M=32, N=16);
plot_surfs(surfs5)
#%%
#add more to chaotic region, most likely these will make it worse.
#these surfs make B^s ~10x larger in chaotic region.
#somewhat expected
rats6 = [(18, 11), (17, 10), (16, 9), (20, 11)]
alist6 = [i[1] for i in rats6]
blist6 = [i[2] for i in rats6]
gl6 = surf_guess(rats5, q_prof)
surfs6 = MID.QFM.new_construct_surfaces(alist6, blist6, gl6, prob, M=32, N=16);
plot_surfs(surfs6)
#%%
curr_surfs = vcat(surfs1, surfs2);
curr_surfs = vcat(surfs1, surfs2, surfs3);
#19/7 surface is overlapping with neighbour.
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4);
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5);

#(11, 7) is cooke d for some reason. Bit surprising as it looks pretty fine.
#just had it twice lol.
curr_surfs = vcat(surfs1, surfs2, surfs3[1:1], surfs3[3:6]);
curr_surfs = vcat(surfs1, surfs2, surfs3[1:1], surfs3[3:6], surfs4);
curr_surfs = vcat(surfs1, surfs2, surfs3[1:1], surfs3[3:6], surfs4, surfs5);
curr_surfs = vcat(surfs1, surfs2, surfs3[1:1], surfs3[3:6], surfs4, surfs5, surfs6);
plot_surfs(curr_surfs, legend=false)

save_object("low_shear_surfs.jld2", curr_surfs)

#%%
rgrid_jac = init_grid(type=:rf, N = 100, start=0.05, stop=0.95)
θgrid_jac = init_grid(type=:af, N = 20) 
ζgrid_jac = init_grid(type=:af, N = 4)
grids_jac = init_grids(rgrid_jac, θgrid_jac, ζgrid_jac)
B, jac, djac = MID.QFM.compute_jac(prob, grids_jac, curr_surfs);
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
ind = find_ind(evals, 0.2373223)
potential_plot(ϕft, grids, ind, label_max=0.05)
#%%

unprob = init_problem(q=q_prof, dens=dens_prof, geo=geo)

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
