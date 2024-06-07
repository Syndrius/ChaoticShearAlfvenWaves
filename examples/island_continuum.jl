

using MID

#testing the island continuum case.

#I guess we should test the no island case first

noA = 1e-3

no_isl = ContIslandT(5, 2, 5/2, 4.0, 0.125, noA)

geo = GeoParamsT(R0=3.0)

#this is a lot of points.
#think something is wrong here tbh.
#unclear wot though
#changing island size did not fix this.
#something wrong with χlist??
#this form is pretty consistent with all of is code.
#theoretically, κ should be within the domain everytime, so what is going on here.
χlistplus = @. -noA - LinRange(0, 1, 801)[2:end]^2 * 0.044
#χlistplus = @. -noA - LinRange(0, 1, 2)[2:end]^2 * 0.044
χlistminus = @. -noA - LinRange(0, 1, 401)[2:end]^2 * 0.00499

#I guess we need to export modedataT..
pmd = MID.ModeDataT(start=0, count=33, incr=2)
tmd = MID.ModeDataT(start=-5, count=3, incr=2)

ω2list_plus = passing_continuum(χlistplus, pmd, tmd, geo, no_isl, 1)
ω2list_minus = passing_continuum(χlistplus, pmd, tmd, geo, no_isl, -1)
#ω2list_plus = compute_continuum_p(geo, isl, χlistplus, 1, 0, -1, 33, 1, 1, 1, compute_toroidal_metric!)|
ψ̄list_plus = ψ̄_p(no_isl, χlistplus, 1)
ψ̄list_minus = ψ̄_p(no_isl, χlistminus, -1)