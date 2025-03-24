
#lets try and see if our interpolation/CT struct is giving the same stuff as Zhisongs code.

using MID

#%%
function qfm_bench_q(r::Float64)

    a = 1.0
    b = 2.0

    return a + b * r^2, 2 * b * r
end


#%%


k = 0.0012

#guess we start in a slab??? maybe that will help the flux surfaces looks uniform?
geo = GeoParamsT(R0=1)
#so the two wave example is given in flux coords, 
#we will just do the equivalent in r coords, may need to change.
#need to actually add islands lol.
#and compare with a poincare plot.
isl1 = IslandT(m0=3, n0=-4, A=k)
isl2 = IslandT(m0=2, n0=-1, A=0.0/2.0)
prob = init_problem(q=qfm_bench_q, geo=geo, met=cylindrical_metric!, isl=isl1, isl2=isl2); 


#find the edge surfaces for constructing the Farey tree.
#bound at (1, 2), r~0.7
#and at (3, 4), r~0.4
plist = [5, 13, 8]
qlist = [8, 21, 13]
slist = 0.5 .* ones(length(plist))
@time surfs = construct_surfaces(plist, qlist, slist, prob);

#looks like good bounding surfaces!
plot_surfs(surfs);

for surf in surfs
    display(surf.ρ)
end

surf_itp = MID.create_surf_itp(surfs);

CT = MID.CoordTsfmT()

MID.coord_transform!(0.62, 2.0, 0.0, CT, surf_itp)

#%%

x = LinRange(0, 10, 100)
y = x .^2

using Plots
plotlyjs()
plot(x, y .^3)
plot(x, y )
