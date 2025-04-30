#what happens when we put surfaces at the boundaries, i.e. will r be greater than 1 and less than 0?
#so we cannot even find the surfaces at the boundaries, we get errors
#this answers the question we need.
#we can just say that the q profile is chosen so that there are two surfaces close to the domain edges.


using MID
using MIDViz
using Plots
using JLD2
#%%

function bounding_q(r::Float64)
    a = 1
    b = 3/2
    return a + b*r^2, 2*b*r
end

R0=4.0

geo = init_geo(R0=R0)

#chaotic case
k = 0.0006
isl = init_island(m0=3, n0=-2, A=k/3)
isl2 = init_island(m0=4, n0=-3, A=k/4)

#prob = init_problem(q=qfm_q, geo=geo, isl=isl, isl2=isl2)#, flr=flr)
prob = init_problem(q=bounding_q, geo=geo, isl=isl, isl2=isl2)#, flr=flr)
#%%
qlist, plist = farey_tree(3, 1, 1, 5, 2)
guess_list = 0.5 .* ones(length(qlist));
#guess_list = @. sqrt(qlist / plist - 0.95) / sqrt(1.5);
#needed for the single island chain case, this one is probably at the sepratrix I guess.
@time surfs = construct_surfaces(plist, qlist, guess_list, prob);
plot_surfs(surfs)

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
jac, djac, jac_tor, djac_tor, coords = MID.QFM.compute_jac()
