
#try to offer some convergence testing for choosing the proper surfaces
#this will include comparison of jacobian and derivatives as more surfaces are added,
#and a comparison of (r, θ, ζ) coord computation (including the interpolation) as more surfaces are added.
#unsure how we are going to quantify this tbh, seems like we would ideally like to have an idea of better or worse jacobain/interpolated values.
#also hard to tell when we are just adding random af surfaces. this is kind of stupid.

#TODO
function rationals(dmax::Int64, qmin::Float64, qmax::Float64)

    #denominator.
    for i in 1:dmax
    end
end

#%%
using MID
using MIDViz
using Plots
using Plots; plotlyjs()
using JLD2
#%%

R0=4.0

geo = init_geo(R0=R0)


#chaotic case
#seems like a good amount for this case as well.
k = 0.0006
isl = init_island(m0=3, n0=-2, A=k/3)
isl2 = init_island(m0=4, n0=-3, A=k/4)
#isl2 = init_island(m0=4, n0=-3, A=0.0)

prob = init_problem(q=qfm_q, geo=geo, isl=isl, isl2=isl2)#, flr=flr)
#%%
#ideally this would be generated form the above function.
qlist = [1, 2, 3, 5, 4, 5, 7, 5, 7, 9, 6, 7, 8, 9, 11, 12, 7, 11, 13];
plist = [1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6];
#think this is a better way of thinking about this.
arr = [1/1, 7/6, 6/5, 5/4, 4/3, 7/5, 3/2, 8/5, 5/3, 7/4, 9/5, 11/6, 2/1, 13/6, 11/5, 9/4, 7/3, 12/5, 5/2]
qarr = (@. sqrt(arr - 239/240) / sqrt(5/3)) #shows that we need more in between 1/1 and 7/6
#also a worry that we have not many surfaces in between the islands.
display(sort(qlist ./ plist))
guess_list = @. sqrt(qlist / plist - 0.995) / sqrt(1.66);
@time surfs1 = construct_surfaces(plist[1:2], qlist[1:2], guess_list[1:2], prob);
@time surfs2 = construct_surfaces(plist[3:4], qlist[3:4], guess_list[3:4], prob);
@time surfs3 = construct_surfaces(plist[5:7], qlist[5:7], guess_list[5:7], prob);
@time surfs4 = construct_surfaces(plist[8:10], qlist[8:10], guess_list[8:10], prob);
@time surfs5 = construct_surfaces(plist[11:16], qlist[11:16], guess_list[11:16], prob);
@time surfs6 = construct_surfaces(plist[17:19], qlist[17:19], guess_list[17:19], prob);
#plot_surfs(surfs)
#%%
scatter(qarr)

#we may want to change our input arrays to have a structure more representive of this, eg array of tuples or something.
arr = [1/1, 26/25, 25/23, 11/10, 7/6, 6/5, 5/4, 4/3, 11/8, 7/5, 10/7, 3/2, 8/5, 5/3, 7/4, 9/5, 11/6, 2/1, 13/6, 11/5, 9/4, 7/3, 12/5, 5/2]
display(arr[1:11:end])
#%%
#these extras seem to give a pretty damn good Jacobain profile.
surfs0 = construct_surfaces([10, 25, 20], [11, 26, 23], [0.2, 0.1, 0.3], prob);
comb_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5, surfs0, surfs6, surfsc);
plot_surfs(comb_surfs)
length(comb_surfs)
#%%
#these ones are in the chaotic region
@time surfsc = construct_surfaces([8, 7], [11, 10], [0.5, 0.5], prob);
#%%
plot_surfs(surfs1)
plot_surfs(surfs2)
plot_surfs(surfs3)
plot_surfs(surfs4)
plot_surfs(surfs5)
plot_surfs(surfs0)
#%%

Nr = 100
Nθ = 50
Nζ = 1
rgrid = init_grid(type=:rf, N = Nr, start = 0.01, stop = 0.99)
θgrid = init_grid(type=:af, N = Nθ, pf=3)
ζgrid = init_grid(type=:af, N = Nζ, pf=-2)

grids = init_grids(rgrid, θgrid, ζgrid)
rgrid_plot, θgrid_plot, _ = MID.Structures.inst_grids(grids);
#%% - #we probably need to start with a more uniform spread of surfaces, and start adding them in between to see what happens.
#this approach seems like an ok start, but we need to have a more sophisticated way of defining the surfaces, or at least a less arbitrary way. At least the farey tree wasn't arbitrary af.
#maybe another option would be to do the faarey tree type thing, but just remove some surfaces if they are considered `to close' to others? Unsure how we would judge that.
total_surfs = [surfs1, surfs2, surfs0, surfsc, surfs3, surfs4, surfs5, surfs6];

curr_surfs = vcat(total_surfs[1], total_surfs[2], total_surfs[3]);
old_jac, old_djac, old_jac_tor, old_djac_tor, old_coords = MID.QFM.compute_jac(prob, grids, curr_surfs);
for i in 4:length(total_surfs)

    curr_surfs = vcat(curr_surfs, total_surfs[i])
    jac, djac, jac_tor, djac_tor, coords = MID.QFM.compute_jac(prob, grids, curr_surfs)

    jac_diff = jac[:, :, 1] .- old_jac[:, :, 1]
    djac_diff = djac[:, :, :, 1] .- old_djac[:, :, :, 1]
    coords_diff = coords[:, :, :, 1] .- old_coords[:, :, :, 1]

    p = contourf(θgrid_plot, rgrid_plot, jac_diff, levels=50, title="Jacobain"*string(i), ylims=(0.1, 0.9))
    display(p)
    p = contourf(θgrid_plot, rgrid_plot, djac_diff[1, :, :], levels=50, title="dJ/dr"*string(i), ylims=(0.1, 0.9))
    display(p)
    p = contourf(θgrid_plot, rgrid_plot, djac_diff[2, :, :], levels=50, title="dJ/dt"*string(i), ylims=(0.1, 0.9))
    display(p)
    p = contourf(θgrid_plot, rgrid_plot, coords_diff[1, :, :], levels=50, title="r"*string(i), ylims=(0.1, 0.9))
    display(p)
    p = contourf(θgrid_plot, rgrid_plot, coords_diff[2, :, :], levels=50, title="t"*string(i), ylims=(0.1, 0.9))
    display(p)
    old_jac = jac
    old_djac = djac
    old_coords = coords
end
#%%


jac, djac, jac_tor, djac_tor, coords = MID.QFM.compute_jac(prob, grids, comb_surfs[1:end-2]);
contourf(θgrid_plot, rgrid_plot, jac[:, :, 1], levels=50)
contourf(θgrid_plot, rgrid_plot, djac[1, :, :, 1], levels=50)
contourf(θgrid_plot, rgrid_plot, djac[2, :, :, 1], levels=50)
contourf(θgrid_plot, rgrid_plot, coords[1, :, :, 1], levels=50)
contourf(θgrid_plot, rgrid_plot, coords[2, :, :, 1], levels=50)




qlist, plist = farey_tree(8, 1, 10, 9, 10)
qlist, plist = farey_tree(3, 4, 5, 9, 10)
display(sort(qlist ./ plist))
length(qlist)
display(maximum(qlist))
display(maximum(plist))
display(sort(qlist))

#from Stuarts/Pers paper, they say they are using farey tree/mediants, but that must just be a guideline, then they are picking smaller rationals to fill gaps/replace larger ones.
arr = [0/1, 1/10, 1/9, 1/8, 2/15, 1/7, 3/20, 2/13, 3/19, 1/6, 3/17, 2/11, 3/16, 1/5, 3/14, 2/9, 3/13, 1/4, 3/11, 5/18, 2/7, 5/17, 3/10, 4/13, 1/3, 4/11, 3/8, 5/13, 2/5, 5/12, 3/7, 7/16, 4/9, 5/11, 1/2, 6/11, 5/9, 9/16, 4/7, 7/12, 3/5, 8/13, 5/8, 7/11, 2/3, 7/10, 12/17, 5/7, 8/11, 3/4, 7/9, 11/14, 4/5, 9/11, 5/6, 11/13, 6/7, 13/15, 7/8, 8/9, 9/10, 1/1]
length(arr)
scatter(sort(qlist ./ plist))

