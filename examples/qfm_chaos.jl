"""
For this one we actually add overlapping islands
Looks like garbage atm, but that is probably just because we dont have enough resolution.
In particular, we don't even have the fourier harmonics of the islands.
Think fff will be better for this.
Shows some indication of working, but needs more surfaces and resolution.
"""

#so investigating the r <0.4 error,
#adding an extra surface at r~0.3, fixes the magnitude of the jacobain, this does not make the spectrum any better.
#considering just teh (3, 2) island, the spectrum is very good, and with the extra surface pretty much perfect.
#just the (4, 3) island, the spectrum is awful, and adding the extra surface does not help
#seems like the (4, 3) island is having some kind of second resonance or something.
#we probbaly need to test a different q-profile to see what is going on tbh.
#really quite confusing what is going on here.
#adding another new surface (15/13) makes a big difference for the second island. still big issues around r=0.25, but r=0.4 seems to have fixed itself.
#these exacmples also show that a single island case may actually be practical.
#adding another surface may have actually fixed it. still unsure about going below 0.2 though!
#but still, this shows that the problem is simple not a dense enough list of surfaces, I guess in this region
#the flux surfaces are changing more than it looks like.
#extra surfaces have fixed the choatic case! not bahd.
#note that the extra surfaces have fixed the chaotic case in the cotninuum
#looks like have also fixed the chaotic case for full spectrum, perhaps needs more investigation.
#but this is very promising! -> need better metric for determing the surfaces and if they are good or not.
#extra surfaces also make the poincare plot in qfm coords work perfectly.
#seems like we needed more surfaces, but we probably still need to be careful around the sepratrix.
#perhaps we need to check that again though! -> things may have changed since we last removed them!

#%%
using MID
using MIDViz
using JLD2

#using Plots; gr()
#using Plots; plotlyjs()

#save_dir = "data/qfm/"

#%%
#define the problem to solve

#should be arbitrary!
R0=4.0

#amp needs further thought!
#define the non-resonant island
#with chaos_q, k=0.0025 is very chaotic while sitll having an inner island chain
#0.0027 seems to be about when there is no structures left
#think k=0.0022 may be the best bet, there are very clear structures but it is also very chaotic elsewhere.
k = 0.0006
isl = init_island(m0=3, n0=-2, A=k/3)
#isl = init_island(m0=3, n0=-2, A=0.0)
isl2 = init_island(m0=4, n0=-3, A=k/4)
#isl2 = init_island(m0=4, n0=-3, A=0.0)
k = 0.008
isl = init_island(m0=3, n0=2, A=k/3)
isl2 = init_island(m0=4, n0=-3, A=0.0)

geo = init_geo(R0=R0)

#to solve non-Hermitian
#flr = MID.Structures.FLRT(δ = 1e-18)

prob = init_problem(q=qfm_q, geo=geo, isl=isl, isl2=isl2)#, flr=flr)
#prob = init_problem(q=qfm_benchmark_q, met=:cylinder, geo=geo, isl=isl, isl2=isl2)#, flr=flr)
un_prob = init_problem(q=qfm_q, geo=geo)


#%%
#Define parameters needed for the poincare plot
Ntraj = 100;
#rlist = collect(LinRange(0.45, 0.65, Ntraj));
rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

#poincare_plot(qfm_benchmark_q, slab_to_plot, Nlaps, Ntraj, 0, 0.0, 0.0, 0.0, R0, isl, MID.Structures.no_isl, rlist)
poincare_plot(prob, Nlaps, Ntraj, rlist,  R0)#, filename=save_dir * "original_poincare.png");

#%%

#testing the straightening of the surfaces at the boundary.

#still need to fix this.
#this does not work, don't think it will ever work tbh. think we can just assume the boundaries at 1.0 and 0.0 are flat.
MM = 4
M = 24
N = 8
bound1 = MID.QFM.straighten_boundary(1.0, MM, M, N, prob)
bound2 = MID.QFM.straighten_boundary(0.000001, MM, M, N, prob)
bound3 = MID.QFM.straighten_boundary(0.6, MM, M, N, prob)


#%%
#Now we construct the qfm surfaces, starting with the 2 boundaries

#unfor that the two cases we want to use don't work...
#so the action function was still set to the test mode...
#this may still work as expected.
##so (1, 1) doesn't work, unsure why, some problem in the fourier transform
#May be that some of the parameters are wrong!
#this is probably good enough, (10, 11) is at ~0.3, (1, 2) is at ~0.8
#@time bounding_surfs = construct_surfaces([1, 1], [2, 3], [0.2, 0.8], prob);
#we may actually be able to use (1, 1) now. doesn't seem like a good idea though, at least without increasing some of the parameters
@time bounding_surfs = construct_surfaces([1, 1], [2, 3], [0.2, 0.8], prob);
#perhaos we shhould be less clever with the bounding surfaces
#(1, 4) is at 0.25 ish
#(3, 4) is at 0.75 ish
#@time bounding_surfs = construct_surfaces([1, 3], [3, 4], [0.5, 0.5], prob);

plot_surfs(bounding_surfs);

#%%
#
#So this is getting stupid af, the args of this need to be in the correct order, otherwise this cooks itself.
#perhaps kwargs are needed for this p/q stuff, or at least make this consistent.
#or perhaps have a check in the construct surfaces function that asserst that q>p or whatever we actually want.
#qlist, plist = farey_tree(4, 11, 10, 3, 1)
#qlist, plist = farey_tree(4, 5, 4, 3, 1)
qlist, plist = farey_tree(4, 1, 1, 2, 1)

#the fifth element is ~1.5, very close to the island, so does not converge.
#deleteat!(qlist, 5)
#deleteat!(plist, 5)
#qlist = [11, 3, 14, 17, 37, 31, 48, 45, 25, 39, 53, 64, 36, 61, 47]
#plist = [10, 1, 11, 12, 25, 23, 35, 34, 21, 32, 43, 53, 31, 52, 41]

#qlist = [1, 4, 3, 3, 2, 5, 5, 5] 
#plist = [1, 3, 1, 2, 1, 2, 3, 4] 

#push!(qlist, 12)
#push!(qlist, 5)
#push!(plist, 11)
#push!(plist, 4)


#naturally, this one will take ages to find.
#this may be to big for this to be practical tbh. At least without improving the alg.
#think the other new surface is worthwhile. Think we can at least generate the spectrum from 0.3 onwards.
#this surface took almost 2 hours to find... not ideal
#this surface has not improved the problem. For some reason
#anything below r=0.4 is completly cooked. 
#adding in this extra surface does not help.
#push!(qlist, 29) #gives a surface for qfm_benchmark at 0.01 ish
#push!(plist, 25)
#push!(qlist, 7) #gives a surface for qfm_benchmark at 0.96 ish
#push!(plist, 2)
#we should probably start doing better guesses, i.e. at least invert the q-profile to compute r for each rational.
#0.99 is because one surface whould lie on 0.0, causing issues.
guess_list = sqrt(0.5) .* sqrt.(qlist ./ plist .- 0.99);
#guess_list[1] = 0.2
#guess_list[2] = 0.8
#deleteat!(qlist, 16)
#deleteat!(plist, 16)
#%%
#For depth of 3, giving 9 surfaces, took ~70s
#fkn stupid af this is plist then qlist, opposite to farey_tree.
#depth of 5 took 33 mins lol. But basically only 2 surfs where the problemo.
#and they are cooked implying a proper solution was not found
@time surfs = construct_surfaces(plist, qlist, guess_list, prob);

#so (1, 1) works now, but goes below r=0. may need to change the q-profile a bit.
plot_surfs(surfs)#, filename=save_dir*"surfs.png");
#%%
extra_surfs = construct_surfaces([15, 13, 20], [16, 15, 21], [0.3, 0.3, 0.26], prob);
comb_surfs = vcat(extra_surfs, surfs);
plot_surfs(comb_surfs)
#%%
#adding the boundary at 0.0 is not good, because the flux surface there is not actually straight. the one at 1.0 is fine though!
#may just need to focus in on the chaotic region.
#ideally we can show the entire plot. In particular, we want to show the eigenfunctions across the whole domain.
push!(surfs, bound1);
push!(surfs, bound2);

#%%

save_object("/Users/matt/phd/MID/data/qfm/chaos_surfs.jld2", surfs)
save_object("/Users/matt/phd/MID/surf_test_isl2.jld2", surfs)

#%%

surfs = load_object("/Users/matt/phd/MID/data/qfm/chaos_surfs.jld2");
#surfs = load_object("surfs5.jld2");
surfs = load_object("/Users/matt/phd/MID/surf_test_isl.jld2");

#%%

#implies we may need a slower growing q-profile, or perhaps it would be better to be able to combine multiple surface runs together, as per Zhisongs add surf.
plot_surfs(surfs)#, filename=save_dir*"surfs.png");


#%%

#now use the qfm surfaces to view the new poincare plot.
##small case with 10, 100 shows pretty stratight B field, this is kinda slow tho!
#this seems to be cooked af now for the chaos case, probably a bad sign for the future.
#think this might actually be broken, or this is more sensitive to the Jacobian being less than ideal.
#so it seems like the real problemo is outside the chaotic region. Pretty annoying tbh.
#perhaps we should stick with the original case, and just consider a truncated domain, 
#that we can finish the draft and then try and fix this garbage later
Ntraj = 40;
#rlist = collect(LinRange(0.4, 0.65, Ntraj));
#so r < 0.4 seems to be completly cooked for every situation...
#not ideal
rlist = collect(LinRange(0.2, 0.8, Ntraj));
Nlaps = 500;

#so this doesn't work v nice.
#guess this tells us that the surfaces are cooked

#this looks to sort of be working now, probably need many more surfaces tbh.
#poincare_plot(qfm_benchmark_q, slab_to_plot, Nlaps, Ntraj, 0, 0.0, 0.0, 0.0, R0, isl, MID.Structures.no_isl, rlist)
x, z = poincare_plot(prob, Nlaps, Ntraj, rlist,  R0, surfs=comb_surfs)#, filename=save_dir * "qfm_poincare.png");

#%%

#need this function, should be inbuilt into MID tbh

function reconstruct_ω(ω, grids)
    #turn a weir doutput thing into a normal output hopefully.
    new_ω = zeros(grids.r.N, grids.θ.N, grids.ζ.N)

    for i in 1:grids.r.N, j in 1:grids.θ.N * grids.ζ.N

        θ, ζ = MID.Indexing.index_to_grid(j, grids)
        new_ω[i, θ, ζ] = ω[i, j]
    end
    return new_ω

end
#%%
#now we construct the cotninuum grids

ϑgrid = MID.Structures.asm_grid(start=0, N=6)#, f_quad=1)
φgrid = MID.Structures.asm_grid(start=-3, N=3)#, f_quad=1)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);

#bounds chosen to not go outside the bounding surfaces
sgrid = MID.Structures.ContGridDataT(N=50, start=0.2, stop=0.8)
cont_grids = init_grids(sgrid, ϑgrid, φgrid)

#%%
#now we construct the continuum
cont_ω = reconstruct_ω(MID.Construct.continuum(prob, cont_grids, surfs), cont_grids);
cont_ω = reconstruct_ω(MID.Construct.continuum(prob, cont_grids, comb_surfs), cont_grids);
#we have removed the perN option, not ideal. should add back in!
#cont_norm_ω = reconstruct_ω(MID.Construct.continuum(un_prob, cont_grids, false), cont_grids);
cont_norm_ω = reconstruct_ω(MID.Construct.continuum(un_prob, cont_grids), cont_grids);


#%%
continuum_plot(cont_ω, cont_grids)#, filename="data/qfm/qfm_cont.png")
continuum_plot(cont_norm_ω, cont_grids)#, filename="data/qfm/unpert_cont.png")



#%%

#now we consider the full spectrum

rgrid = MID.Structures.rfem_grid(N=80, start=0.2, stop=0.8)
θgrid = MID.Structures.asm_grid(start=0, N=6, f_quad=4)
ζgrid = MID.Structures.asm_grid(start=-3, N=3, f_quad=4)#, incr=2)
#rgrid = MID.Structures.rfem_grid(N=50, start=0.3, stop=0.8)
#θgrid = MID.Structures.afem_grid(N=8, pf=3)
#ζgrid = MID.Structures.afem_grid(N=2, pf=-2)#, incr=2)

grids = init_grids(rgrid, θgrid, ζgrid)

#%%

#solver = init_solver(nev=100, target = 0.3, prob=prob)
solver = init_solver(nev=100, targets = [0.1, 0.2, 0.3, 0.4], prob=prob)

#%%

#seems to be a fair bit slower!
#makes sense, will want to use MIDParallel oneday.
evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, solver=solver, surfs=comb_surfs);

evals_norm, ϕ_norm, ϕft_norm = compute_spectrum(prob=un_prob, grids=grids, full_spectrum=false);

evals_weird, ϕ_weird, ϕft_weird = compute_spectrum(prob=prob, grids=grids, full_spectrum=false);

#%%

continuum_plot(evals)#, filename="data/qfm/qfm_spectrum.png")
continuum_plot(evals_norm)#, filename="data/qfm/norm_spectrum.png")
continuum_plot(evals_weird)#, filename="data/qfm/weird_spectrum.png")

#%%

ind = find_ind(evals, 0.28715)
ind = 10

#hard to be certain, but looks a lot better now, may even be working!
potential_plot(ϕft, grids, ind)
