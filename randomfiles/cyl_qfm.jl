
#need to decide on a q-profile to use for the cylindrical case for qfm.
#perhaps we can now do the stupid RMP high shear case
#although that will most likely still be fkn awful

#we probably instead want to do a q-profile that varies between 1 and 2 or something
#given we no longer care about global modes (RIP!) low shear may not be as important, however we still don't want shear to be to high because then the continuum becomes cooked af.
#note that we may actually want to move to a slab
#problem then will be the annoying normalisation of the eigenvalues
#but that may be easier to fix than dealing with the axis.
#axis does completly cook the flux surfaces
#does seem like if the flux surfaces for r<chaos where more like the r>chaos surfaces we would have less problemos
#slab doesn't seem to fix this
#perhaps we need to work with flux coordinates
#or juust remove the r factor in Bθ?
#think removing r is the same as working in ψ
#think that should probably be our next go.
#this will obviously be a complete fkn nightmare with qfm
#but we are really grasping at straws.
using MID
using MIDViz
using Plots; plotlyjs()
using Plots; gr()
#%%

function cyl_qfm_q(r::Float64)

    return 1.0 + r^2, 2*r
end
#%%
geo = init_geo(R0=10.0)

prob = init_problem(q = cyl_qfm_q, geo=geo, met=:slab)

#%%
rgrid = init_grid(type=:rc, N=100)
θgrid = init_grid(type=:as, N=5, start=1)
ζgrid = init_grid(type=:as, N=2, start=-2)

grids = init_grids(rgrid, θgrid, ζgrid)

#%%

evals = compute_continuum(prob, grids);

continuum_plot(evals, grids)
#%%
#poincare plot

isl1 = init_island(m0=3, n0=-2, A=0.000)
isl2 = init_island(m0=7, n0=-5, A=0.003)
isl3 =  init_island(m0=8, n0=-5, A=0.002)
isl4 =  init_island(m0=6, n0=-5, A=0.001)

isls = [isl1, isl2, isl3, isl4]
prob = init_problem(q = cyl_qfm_q, geo=geo, met=:slab, isls=isls)
#%%

Ntraj = 100;
#rlist = collect(LinRange(0.45, 0.65, Ntraj));
rlist = collect(LinRange(0.01, 1.0, Ntraj));
Nlaps = 500;

#poincare_plot(qfm_benchmark_q, slab_to_plot, Nlaps, Ntraj, 0, 0.0, 0.0, 0.0, R0, isl, MID.Structures.no_isl, rlist)
poincare_plot(prob, Nlaps, Ntraj, rlist)#, filename=save_dir * "original_poincare.png");
