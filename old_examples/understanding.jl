
#having another go at understanding TAEs, and see if we can pinpoint island gap modes/continuum modes
#and see if we can see any modification to TAEs because of the island


#the two types of TAEs we are seeing seem to reflect the top and bottom of the gap nicely.
#bottom of the gap matches our usual tae mode strutcure, where m=2/3 modes are added together, ie both mode contributions have the same sign
#whereas the top case modes are of opposite signs.

#still unclear why this would be...
#seems like we have narrowed the problem to this.

#part of the story is given in Cheng's paper
#the top gaps case, where the modes are subtracted causes the net wave to peak at θ=π, ie the inner side of the torus where the magnetic field is stronger, which may explain the higher frequency.
#similary, we get the opposite for the bottom gap case, where the net wave is larger at θ=0,2π where the magnetic field is weaker.

#this is not the full story but is another peice of the ol puzzle.


#seems like there are two possible ways the waves can be combined at the gap location. They can be added or subtracted. Adding them results in a lower frequency wave, due to the amplitude being higher when B is lower, or the waves can be subtracted, resulting in a higher frequency wave as the amplitude is larger when B is larger.

#this can be viewed as the contrasting ways the waves can interfer, i.e. in phase or out of phase, which can only occur once toroidaicity allows the two waves to interact.
#this kind of ruins our previous argument though, because we justified toroidicity modification by requiring that the potential be modified according to the poloidal variation. Perhaps this aspect is more complicated?
#i.e perhaps it is not so simple that the potential should increase where B is larger?
#perhaps it is just that we expect the potential to be modified at θ=0,π,2π? perhaps the potential could equally decrease? This could be similar to our previous args, in that the periodic perturbation could chaneg the direction of B a bit, which would actually decrease k⋅B even if mag B is getting bigger?
#ideally we can show this, because that might be enough to craft the full story.
#think we need to return to our understanding of why frequency decreases ever!

#feel like we can explain the modification of normal modes as the `requirement` for modification to adjust to the new feild.
#But the gap parts need to be explained as if both waves already exist and interfere with each other
#hard to combine both of these explainations.

#ie does toroidicty perturbtaion `allow` for coupling (suddenly!) maybe more importantly,
#why don't they interfer with each other before.
#or are modes modified by neighbours to correct for the torodicty
#I think it is both but also neither.


#perhaps lower frequency modes need less of a correction
#this kinda makes sense, as they would maybe see `less` of the perturbation??
#this does seem to be pretty consistent so may need to undertsand this...
#seems like the proportion of change in B would still be the same though??

#perhaps in cylindrical case, the normal modes are true fourier harmonics
#but now the normal modes are made up of multiple harmonics -> perhaps this is what now allows them to interact??? i.e individual fourier harmonics are indep (net overlap is zero!) 
#but when we have two modes at similar freqs (i.e gaps!) that are capable of interfering,
#we find two new normal modes. i.e the top and the bottom!
#Still don't really think this explains why some evals are sums while some are subtracts.
#looks like this is because any time a mode is modified it is because it is appraoching a gap!
#or it is part of a gap, even if it is a long way away from the gap part
#consider m=1,2 and n=-2 in cylindrical and toroidal cases, we still see a widening of the distance between them!

#not sure this is the full story but we are getting closer!
#question still remaining. Just because we can add the fourier modes, why would this happen? Why would the potential get larger where B is smaller? That seems backwards? i.e how can the electric potential be largest in the region where B is smallest? i.e How can this wave exist? Seems like it must be travelling against B or something??
#chief difference seems to be sign of ϕ at θ=π (max of B) If sign of ϕ matches sign of B then |ϕ| decreases, else |ϕ| increases...

#so one solution is like ϕ + B the other is like ϕ - B. (but with relative magnitude of B.)
#I think the two options are co- and counter-propogating. We will need to show this properly!
#Bill's review states quite clearly that they are counter propogating. Need to understand their propogating nature,
#this may also help us undertsand tae's as they are standing waves? right...?

#I guess gaps can only ever occur if one is going backwards and one forwards, otherwise both modes would have the same behaviour with r, ie both be increasing or decreasing with r.

#is this perhap all the peices?? maybe... start trying to write this up and we will almost certainly find that it is not.

function test_q(r)
    #test case for rsae from Bill's review
    q = 1 - 0.5 * r + r^3
    dq = -0.5 + 3 * r^2
    return q, dq
end


using MID

using Plots; plotlyjs()

N = 300;
grids = init_grids(N=N, mstart=-10, mcount=26, nstart=-10, ncount=1, nincr=4);

isl = IslandT(A=0e-5, m0=5, n0=4);
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=island_damping_q, isl=isl, geo=geo);
#prob = init_problem(q=fu_dam_q, isl=isl, geo=geo);
#prob = init_problem(q=test_q, isl=isl, geo=geo);


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true, σ=(0.395/geo.R0)^2, nev=100);


reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)#, ymax=10)

#0.376 for bottom of gap
#0.423 for top of gap -> with N=200, I think that is another TAE thingo.
#this is a separate issue! -> it appears in the other code as well!, it is not as clearly in the gap though!
#try 0.43 for top
tae_ind = find_ind(ω, 0.183)
tae_ind = find_ind(ω, 0.66)

tae_ind = find_ind(ω, 0.1518)
tae_ind = find_ind(ω, 0.383) #why is this one basically a cylindrical mode??? this occurs at large r so should have toroidal corrections?

tae_ind = find_ind(ω, 3.576557)
tae_ind = find_ind(ω, 4.46)

display(ω[tae_ind])

#for R0=3, wtf is the gap mode at ω≈0.44, is it a tae? we have seen similar stuff before (with islands)
#this is when we get an antisymmetric pair of m=2 and m=3.
plot_potential(grids=grids, ϕ=ϕ, ind=tae_ind, n=1)


z = construct_surface(ϕ, length(ω), grids, 1);

plot_surface(z, grids, tae_ind)



ω_cont = continuum(prob=prob, grids=grids);
#this doesn't tell me which n is which fk sake.
plot_continuum(ω = ω_cont, grids=grids, n=-6, ymax=1)#, filename="data/continuum.png")
#in cylindrical limit,

#ω1^2 = (m/q + n)^2 with some R's or someshit.
#ω2^2 = ((m+1)/q + n)^2

#q(0.5) = 1.25

#so ω1^2 = (2/1.25 - 2)^2 -> ω1 = -0.4
#ω2^2 = (3/1.25 - 2)^2 -> ω2 = 0.4

#two beats are \abs(f1-f2)/2

#doesn't make any sense, beats should have the same amplitude and same frequency. think this is wrong approach.

island_damping_q(0.5)
