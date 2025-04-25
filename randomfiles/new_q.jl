#new q-profile for qfm, chosen so that we can have qfm surfaces very close to the boundary, allowing us to compute the spectrum across the full domain.

using MID
using MIDViz
using Plots
using JLD2
#%%
function new_q(r::Float64)
    a = 0.995833333
    b = 1.66666666
    return a + b * r^2, 2 * b * r
end

R0=4.0

geo = init_geo(R0=R0)


#chaotic case
#seems like a good amount for this case as well.
k = 0.0006
isl = init_island(m0=3, n0=-2, A=k/3)
isl2 = init_island(m0=4, n0=-3, A=k/4)
#isl2 = init_island(m0=4, n0=-3, A=0.0)

prob = init_problem(q=new_q, geo=geo, isl=isl, isl2=isl2)#, flr=flr)

#%%
Ntraj = 100
rlist = collect(LinRange(0.0001, 1.0, Ntraj));
rlist = collect(LinRange(0.35, 0.65, Ntraj));
Nlaps = 500
poincare_plot(prob, Nlaps, Ntraj, rlist, R0);
#%%

#seems like a good concentration of surfaces for r>0.5
#but nowhere near enough for below r=0.5
#starting at 1/1 and 5/2 seems bad, missing way to much of the low r values.
#starting at 1/1 and 2/1 is better, but this is still kind of a shit way to do this.
qlist, plist = farey_tree(4, 1, 1, 2, 1)
#probably want a function for this, just iterates the denominator and takes all numerators that are between 0.0 and 2.5 (edges of q-profile) and makes sure it is not already in the array.
qlist = [1, 2, 3, 5, 4, 5, 7, 5, 7, 9, 6, 7, 8, 9, 11, 12]
plist = [1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5]
guess_list = @. sqrt(qlist / plist - 0.995) / sqrt(1.66);
@time surfs = construct_surfaces(plist, qlist, guess_list, prob);
plot_surfs(surfs)

#%%
#these extras seem to give a pretty damn good Jacobain profile.
extra_surfs = construct_surfaces([10, 25, 20], [11, 26, 23], [0.2, 0.1, 0.3], prob);
comb_surfs = vcat(surfs, extra_surfs);
plot_surfs(comb_surfs)

#%%
#now the qfm poincare plot
Ntraj = 40
#starting this below 0.4 causes immediate issues for the Jacobain being enourmous/negative
rlist = collect(LinRange(0.4, 0.9, Ntraj));
Nlaps = 500
x, z = poincare_plot(prob, Nlaps, Ntraj, rlist, R0, surfs=surfs_chaos);

#%%
Nr = 100
Nθ = 20
Nζ = 1
rgrid = init_grid(type=:rf, N = Nr, start = 0.01, stop =0.99)
θgrid = init_grid(type=:af, N = Nθ, pf=5)
ζgrid = init_grid(type=:af, N = Nζ, pf=-2)

grids = init_grids(rgrid, θgrid, ζgrid)
jac, djac, jac_tor, djac_tor = MID.QFM.compute_jac(prob, grids, surfs);
jac, djac, jac_tor, djac_tor = MID.QFM.compute_jac(prob, grids, comb_surfs);
#%%
rgrid_plot, θgrid_plot, _ = MID.Structures.inst_grids(grids);
contourf(θgrid_plot, rgrid_plot, jac[:, :, 1], levels=50)
contourf(θgrid_plot, rgrid_plot, djac[1, :, :, 1], levels=50)
contourf(θgrid_plot, rgrid_plot, djac[2, :, :, 1], levels=50)
