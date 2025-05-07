#testing newest surface construction

using MID
using MIDViz
using Plots; plotlyjs()


#%%

function q_prof(r::Float64)
    #we are going to stick with the q-prof
    #n=-3 harmonics for perturbations
    #if we need to look at tae, we will need to change the density profile
    
    a = 1.05
    b = 0.5
    c = 9.0
    d = 0.0
    return a + b*r^2 + c*r^6, 2 * b * r + 6*c*r^5
    #return a + b*r^2 + c*r^10 + d*r^4, 2*b*r + 10*c*r^9 + 4*d*r^3
    #return a + b*r^6, 6 * b * r^5
    #return a + b*r^4, 4 * b * r^3
end
function dens_prof(r::Float64)

    return 1.0 #- tanh((r-0.85)/0.15)
end
#%%

rgrid = init_grid(type=:rf, N=60)
θgrid = init_grid(type=:as, start=0, N=9)
ζgrid = init_grid(type=:as, start=-5, N=5)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%

solver = init_solver(nev = 200, targets=[0.25, 0.30, 0.35], prob=prob)

#%%
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);
#%%

k = 0.0015
isl1 = init_island(m0=7, n0=-3, A=k/7)
isl2 = init_island(m0=8, n0=-3, A=k/8)
isls = [isl1, isl2]

geo = init_geo(R0=4.0)
prob = init_problem(q=q_prof, dens=dens_prof, geo=geo, isls=isls)
#%%
#last one for depth 7 didn't seem to work, 16/7 I think.
rats = lowest_rationals(6, q_prof(0.0)[1], q_prof(0.7)[1])
alist = [i[1] for i in rats]
blist = [i[2] for i in rats]
gl = 0.5 .* ones(length(alist))
#%%
@time og_surfs = construct_surfaces(blist[1:end], alist[1:end], gl[1:end], prob);

#so much faster, depth 6 is over 2mins vs 15s.
#still needs to be fixed, and there may be some small improvements to be made.
@time new_surfs = MID.QFM.new_construct_surfaces(alist[1:end], blist[1:end], gl[1:end], prob);

#%%

plot_surfs(og_surfs)
plot_surfs(new_surfs)

#%%

arr1 = zeros(7, 4)
arr2 = zeros(7, 4) .+ 1
arr3 = zeros(7, 1) .+ 4

[arr1  arr2  arr3]

