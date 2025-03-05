
using MID
using JLD2

#using Interpolations

#so for this we probably just want the interpolant right?
#we don't actuallly care about the `true` surface object?
#surely we can use the interpoaltion to plot if needed?
#perhaps we start separate, but then combine if it is pointless to have two objects.
using Plots; plotlyjs()

function two_wave_q(r::Float64)

    #flux q is 1/r.
    q = 2/ r^2
    dq = -4 /r^3
    return q, dq
end


function qfm_q(r::Float64)

    #chosen so 2/1 isl at 0.65 and 3/2 isl at 0.5
    a = 0.7754
    b = 2.8986
    return a + b*r^2, 2 * b * r

end


k = -0.0018

geo = GeoParamsT(R0=1)
#so the two wave example is given in flux coords, 
#we will just do the equivalent in r coords, may need to change.
#need to actually add islands lol.
#and compare with a poincare plot.
isl1 = IslandT(m0=3, n0=-2, A=k/3.0)
isl2 = IslandT(m0=2, n0=-1, A=k/2.0)
prob = init_problem(q=two_wave_q, geo=geo, met=slab_metric!, isl=isl1, isl2=isl2); 

#prob = init_problem(q=two_wave_q, geo=geo, met=cylindrical_metric!, isl=isl1, isl2=isl2); 

#in theory this won't need the actual grids used to solve
#but it will need some kind of `grid` type structure that tells how big the surface approximation should be etc.
#looks like this is working, just slow af.
#wthis will need a serious zooming, especially if we require computing more than ~5 surfaces
#or if we add the action gradient to this.
#this surfs obj is quite large as well, so ideally this would not be stored long term.
#rather, just used to create the interpolants
#surfs = construct_surfaces([5, 13, 8, 11, 3], [8, 21, 13, 18, 5], prob);

8/13
13/21
5/8

#think these rationals will be different for flux vs r.
#this is getting out of hand.
#guess we can see why Zhisong had a save to file function
#probably worth implementing.
#should be easy enough with jld2?
surfs = construct_surfaces([8, 13, 5], [13, 21, 8], prob);

save_object("qfmTest/test_surfs.jld2", surfs)




#this is a wee bit different to before, but the magnetic field is just a wee bit different so it is to be expected.
#although this is really quite different...
#less peaks, not as flat etc
#ok so we may just need to come up with our own example case?
#this flux one is not working.
#ok we are back!
MID.QFM.plot_surfs(surfs2);


surfs = load_object("qfmTest/test_surfs.jld2")
#this at least seems to be very quick despite the stupidity.
surf_itp = MID.QFM.create_surf_itp(surfs);

CT = MID.QFM.CoordTsfmT()

MID.QFM.coord_transform(1.1095, 0.2, 0.0, CT, surf_itp)


#prob = init_problem(q=test_q, geo=geo)#, met=no_delta_metric!); 

#rgrid = init_fem_grid(N=Nr)
θgrid = asm_grid(start=-1, N=5)
ζgrid = asm_grid(start=-2, N=5)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);

rgrid = MID.ContGridDataT(N=10, start=1.11, stop=1.115)
grids = init_grids(rgrid, θgrid, ζgrid)

#ok so working, in that there are no errors
#all returned evals are wopping though, ~10^20
#this could be a numerical problem or a problem with the 
#example case, probbaly both
#we will need to fix the interpolation outside the flux surfaces
#probably easiest to just add a flux surface at r=0 and r=1 tbh.
ω = qfm_continuum(prob, grids, surfs)

