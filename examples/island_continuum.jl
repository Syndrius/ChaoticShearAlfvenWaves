

using MID
using Plots

#testing the island continuum case.

#I guess we should test the no island case first

noA = 1e-10

#this is being constructed in the wrong order!!! Fixes domain issue.
no_isl = ContIslandT(5, 2, noA, 5/2, 4.0, 0.125)

geo = GeoParamsT(R0=3.0)

#this is a lot of points.
#think something is wrong here tbh.
#unclear wot though
#changing island size did not fix this.
#something wrong with χlist??
#this form is pretty consistent with all of is code.
#theoretically, κ should be within the domain everytime, so what is going on here.

#this is far to many points for us to be using tbh.
#χlistplus = @. -noA - LinRange(0, 1, 801)[2:end]^2 * 0.044
χlistplus = @. -noA - LinRange(0, 1, 401)[2:end]^2 * 0.044
#χlistplus = @. -noA - LinRange(0, 1, 2)[2:end]^2 * 0.044
χlistminus = @. -noA - LinRange(0, 1, 401)[2:end]^2 * 0.00499


#display(@. (-χlistplus[1:5]+noA) / (2*noA))

#I guess we need to export modedataT..
pmd = MID.ModeDataT(start=0, count=33, incr=1)
tmd = MID.ModeDataT(start=-5, count=3, incr=2)

#ω2list_plus = passing_continuum(χlistplus, pmd, tmd, geo, no_isl, 1);
ω2list_plus = island_continuum(χlistplus, pmd, tmd, geo, no_isl, 1);
#ω2list_minus = passing_continuum(χlistminus, pmd, tmd, geo, no_isl, -1);
ω2list_minus = island_continuum(χlistminus, pmd, tmd, geo, no_isl, -1);
#ω2list_plus = compute_continuum_p(geo, isl, χlistplus, 1, 0, -1, 33, 1, 1, 1, compute_toroidal_metric!)|
ψ̄list_plus = MID.IslandContinuum.compute_ψ̄(no_isl, χlistplus, 1)
ψ̄list_minus = MID.IslandContinuum.compute_ψ̄(no_isl, χlistminus, -1)


#33 * 1 is mcount by ncount, should be done more automatically in the future
rp = repeat(sqrt.(2 .* ψ̄list_plus), 1,  pmd.count * tmd.count)
rm = repeat(sqrt.(2 .* ψ̄list_minus), 1, pmd.count * tmd.count)

#ayyy we have the very basic funcitonality back!!!!
#now we just need to understand how it works and use it properly... Simple!
#plotting this with some proper colours would be nice.
scatter(rp, sqrt.(abs.(ω2list_plus .* geo.R0^2)), ylimits=(0, 1), legend=false)
scatter!(rm, sqrt.(abs.(ω2list_minus .* geo.R0^2)))




#now small island case

smallA = 1e-5

small_isl = ContIslandT(5, 2, smallA, 5/2, 4.0, 0.125)

geo = GeoParamsT(R0=3.0)

#this is far to many points for us to be using tbh.
#srs wot the fk are these ranges???
#seem completly random..
χlistplus = @. -smallA - LinRange(0, 1, 401)[2:end]^2 * 0.005
#χlistplus = @. -noA - LinRange(0, 1, 2)[2:end]^2 * 0.044
χlistminus = @. -smallA - LinRange(0, 1, 401)[2:end]^2 * 0.003


pmd = MID.ModeDataT(start=-20, count=41, incr=1)
tmd = MID.ModeDataT(start=-5, count=6, incr=2)

#ω2list_plus = passing_continuum(χlistplus, pmd, tmd, geo, small_isl, 1);
ω2list_plus = island_continuum(χlistplus, pmd, tmd, geo, small_isl, 1);
#ω2list_minus = passing_continuum(χlistminus, pmd, tmd, geo, small_isl, -1);
ω2list_minus = island_continuum(χlistminus, pmd, tmd, geo, small_isl, -1);
#ω2list_plus = compute_continuum_p(geo, isl, χlistplus, 1, 0, -1, 33, 1, 1, 1, compute_toroidal_metric!)|
ψ̄list_plus = MID.IslandContinuum.compute_ψ̄(small_isl, χlistplus, 1)
ψ̄list_minus = MID.IslandContinuum.compute_ψ̄(small_isl, χlistminus, -1)

rp = repeat(sqrt.(2 .* ψ̄list_plus), 1,  pmd.count * tmd.count)
rm = repeat(sqrt.(2 .* ψ̄list_minus), 1, pmd.count * tmd.count)

#ayyy we have the very basic funcitonality back!!!!
#now we just need to understand how it works and use it properly... Simple!
#plotting this with some proper colours would be nice.
scatter(rp, sqrt.(abs.(ω2list_plus .* 3^2)), ylimits=(0, 0.35), legend=false, xlimits=(0.3, 0.7))
scatter!(rm, sqrt.(abs.(ω2list_minus .* 3^2)))



#`normal` island case


A = 1e-4

isl = ContIslandT(5, 2, A, 5/2, 4.0, 0.125)

geo = GeoParamsT(R0=3.0)

#this is far to many points for us to be using tbh.
#srs wot the fk are these ranges???
#seem completly random..
χlistplus = @. -A - LinRange(0, 1, 401)[2:end]^2 * 0.005
#χlistplus = @. -noA - LinRange(0, 1, 2)[2:end]^2 * 0.044
χlistminus = @. -A - LinRange(0, 1, 401)[2:end]^2 * 0.003


pmd = MID.ModeDataT(start=-20, count=41, incr=1)
tmd = MID.ModeDataT(start=-5, count=6, incr=2)

#ω2list_plus = passing_continuum(χlistplus, pmd, tmd, geo, isl, 1);
ω2list_plus = island_continuum(χlistplus, pmd, tmd, geo, isl, 1);
#ω2list_minus = passing_continuum(χlistminus, pmd, tmd, geo, isl, -1);
ω2list_minus = island_continuum(χlistminus, pmd, tmd, geo, isl, -1);
#ω2list_plus = compute_continuum_p(geo, isl, χlistplus, 1, 0, -1, 33, 1, 1, 1, compute_toroidal_metric!)|
ψ̄list_plus = MID.IslandContinuum.compute_ψ̄(isl, χlistplus, 1)
ψ̄list_minus = MID.IslandContinuum.compute_ψ̄(isl, χlistminus, -1)

rp = repeat(sqrt.(2 .* ψ̄list_plus), 1,  pmd.count * tmd.count)
rm = repeat(sqrt.(2 .* ψ̄list_minus), 1, pmd.count * tmd.count)

#ayyy we have the very basic funcitonality back!!!!
#now we just need to understand how it works and use it properly... Simple!
#plotting this with some proper colours would be nice.
scatter(rp, sqrt.(abs.(ω2list_plus .* 3^2)), ylimits=(0, 0.35), legend=false, xlimits=(0.3, 0.7))
scatter!(rm, sqrt.(abs.(ω2list_minus .* 3^2)))


#now we gotta do the inside of the island.
#this is what Zhisong does, takes a while though!
#pmd = MID.ModeDataT(start=-50, count=101, incr=1)
#tmd = MID.ModeDataT(start=-11, count=12, incr=2)

#now consider less modes to see what happens.
#with lower mode numbers, the gaps still seems to be the same
#good for us, hopefully means we can test gap width more quickly.
#kind of seems like the behaviour further from the gap gets more and more cooked.
pmd = MID.ModeDataT(start=-30, count=61, incr=1)
tmd = MID.ModeDataT(start=-8, count=9, incr=2)

#@. doesn't work with exp of linrange for some reason!
χlist1 =  1 .- 0.01 * exp.( -1 .* collect(LinRange(0, 12, 11)))
χlist2 = LinRange(0, sqrt(χlist1[1]), 191)[1:end-1] .^2 .+ 1e-5
#ok this matches Zhisong, again no fkn idea what the hell this list is.
χlist = A .- vcat(χlist2, χlist1) .* 2 .* A

#may actually be different!
#not ideal, will need to test with less modes. takes like ~5mins.
#ω2list = trapped_continuum(χlist, pmd, tmd, geo, isl);
ω2list = island_continuum(χlist, pmd, tmd, geo, isl, 0);

#display(ω2list)


width = 4 * sqrt(isl.A * isl.q0^2/isl.qp)
ψ_isl = 2 * width / (π * isl.m0)
ψ̄m = ψ_isl

ψ̄list = MID.IslandContinuum.compute_ψ̄(isl, χlist, 0)
r = repeat(sqrt.(2 .* ψ̄list), 1,  pmd.count * tmd.count)


#eyy looks like Zhisongs, pretty unclear what any of it means tbh!
scatter(r, sqrt.(abs.(ω2list .* geo.R0^2)), legend=false, markersize=0.1, ylimits=(0, 0.3))



#I guess we should now consider our island...
#will need to understand wot the hek is going on with χ lol
#may also want R0 to be smaller so gap is larger

#ie 1e-5 is the part we pass into our other code.
#so with 1e-6 we actually expect the tae to not experience damping
#but 1e-5 perhaps not? there is like 1 `strand` that it could interact with. But majority is below.
#similar with 5e-5, perhaps we just don't have enough modes?
A = 1e-4 / 4 #4 to reflect r(1-r) factor

#so with small islands the structures are almost independent?

#educated guess lol.
isl = ContIslandT(5, 4, A, 5/4, 0.8, 0.125)

geo = GeoParamsT(R0=10.0)

#start with this???
pmd = MID.ModeDataT(start=-30, count=61, incr=1)
tmd = MID.ModeDataT(start=-18, count=8, incr=4)
#tmd = MID.ModeDataT(start=-8, count=10, incr=2)

#χlist = LinRange(-A+A*0.05, A-A*0.05, 50) #this is fked, think we need to cluster near the spratrix.

#just try what Zhisong did?
χlist1 =  1 .- 0.01 * exp.( -1 .* collect(LinRange(0, 12, 11)))
χlist2 = LinRange(0, sqrt(χlist1[1]), 191)[1:end-1] .^2 .+ 1e-5
#ok this matches Zhisong, again no fkn idea what the hell this list is.
χlist = A .- vcat(χlist2, χlist1) .* 2 .* A

ω2list = trapped_continuum(χlist, pmd, tmd, geo, isl);

#width = 4 * sqrt(isl.A * isl.q0^2/isl.qp)
#ψ_isl = 2 * width / (π * isl.m0)
#ψ̄m = ψ_isl
ψ̄list = MID.IslandContinuum.ψ̄_t(isl, χlist)
r = repeat(sqrt.(2 .* ψ̄list), 1,  pmd.count * tmd.count)

#we should normalise the x-axis, it is going from centre of island to edge I think.
scatter(r, sqrt.(abs.(ω2list .* geo.R0^2)), legend=false, markersize=0.1, ylimits=(0.3, 0.5))#.2, 0.6))

hline!([0.425, 0.425]) #top, 
hline!([0.3764, 0.3764]) #bottom
hline!([0.380, 0.380]) #tae

#our gap has max at 0.424
#min at 0.3764
#tae at ~0.380




A = 1e-4

isl = ContIslandT(5, 2, A, 5/2, 4.0, 0.125)

geo = GeoParamsT(R0=10.0)


pmd = MID.ModeDataT(start=-10, count=2, incr=1)
tmd = MID.ModeDataT(start=-6, count=2, incr=2)

#@. doesn't work with exp of linrange for some reason!
χlist1 =  1 .- 0.01 * exp.( -1 .* collect(LinRange(0, 12, 2)))[1:end-1]
χlist2 = LinRange(0, sqrt(χlist1[1]), 2)[1:end-1] .^2 .+ 1e-5
#ok this matches Zhisong, again no fkn idea what the hell this list is.
χlist = A .- vcat(χlist2, χlist1) .* 2 .* A

#may actually be different!
#not ideal, will need to test with less modes. takes like ~5mins.
ω2list = island_continuum(χlist, pmd, tmd, geo, isl, 0);

ω2list_old = passing_continuum(χlist, pmd, tmd, geo, isl);

#seeing tiny differences, I think form actually computing gu vs taking inverse.
display(ω2list)# .- 
display(ω2list_old)



χlistplus = @. -A - LinRange(0, 1, 401)[2:end]^2 * 0.005
#χlistplus = @. -noA - LinRange(0, 1, 2)[2:end]^2 * 0.044
χlistminus = @. -A - LinRange(0, 1, 401)[2:end]^2 * 0.003




#I guess we need to export modedataT..
pmd = MID.ModeDataT(start=0, count=20, incr=1)
tmd = MID.ModeDataT(start=-5, count=10, incr=2)

ω2list_plus_old = passing_continuum(χlistplus, pmd, tmd, geo, isl, 1);
#ω2list_minus_old = passing_continuum(χlistminus, pmd, tmd, geo, isl, -1);


ω2list_plus = island_continuum(χlistplus, pmd, tmd, geo, isl, 1);
#ω2list_minus = island_continuum(χlistminus, pmd, tmd, geo, isl, -1);


#thinkt they are the same now, hard to be certain though!
display(ω2list_plus .- ω2list_plus_old)
