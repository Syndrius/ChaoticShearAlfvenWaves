"""

This file shows a simple example of the qfm surfaces. Here we add a non-resonant perturbation (in the same form as our magentic islands) and compare the spectrum computed with and without qfm surfaces and also compare these to the unperturbed spectrum.
"""
#%%
using Revise
using MID
using MIDViz
using JLD2

#using Plots; gr()
using Plots; plotlyjs()


#%%
#define the problem to solve

R0=4.0

#amp needs further thought!
#define the non-resonant island
isl = init_island(m0=3, n0=2, A=0.005)

geo = init_geo(R0=R0)

prob = init_problem(q=qfm_benchmark_q, geo=geo, isl=isl)
un_prob = init_problem(q=qfm_benchmark_q, geo=geo)


#%%
#Define parameters needed for the poincare plot
Ntraj = 80;
rlist = collect(LinRange(0.25, 0.75, Ntraj));
Nlaps = 500;

#poincare_plot(qfm_benchmark_q, slab_to_plot, Nlaps, Ntraj, 0, 0.0, 0.0, 0.0, R0, isl, MID.Structures.no_isl, rlist)
poincare_plot(prob, Nlaps, Ntraj, rlist,  R0, filename="data/qfm/bench_og_poincare.png");

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
qlist, plist = farey_tree(5, 2, 1, 3, 1)


guess_list = 0.5 .* ones(length(qlist));
guess_list[1] = 0.2
guess_list[2] = 0.8

#%%
#For depth of 3, giving 9 surfaces, took ~70s
#fkn stupid af this is plist then qlist, opposite to farey_tree.
@time surfs = construct_surfaces(plist, qlist, guess_list, prob);

#%%

save_object("data/qfm/benchmark_surfs.jld2", surfs)

#%%

surfs = load_object("data/qfm/benchmark_surfs.jld2");

#%%

#implies we may need a slower growing q-profile, or perhaps it would be better to be able to combine multiple surface runs together, as per Zhisongs add surf.
#still a big gap at the edge unfort, perhaps in general we define the domain based on where the surfaces are? Rather than trying to put the surfaces where we want?
plot_surfs(surfs, filename="data/qfm/benchmark_surfs.png");

#%%

#now use the qfm surfaces to view the new poincare plot.
##small case with 10, 100 shows pretty stratight B field, this is kinda slow tho!
Ntraj = 40;
rlist = collect(LinRange(0.25, 0.75, Ntraj));
Nlaps = 300;

#poincare_plot(qfm_benchmark_q, slab_to_plot, Nlaps, Ntraj, 0, 0.0, 0.0, 0.0, R0, isl, MID.Structures.no_isl, rlist)
x, z = poincare_plot(prob, Nlaps, Ntraj, rlist,  R0, surfs=surfs, filename="data/qfm/qfm_poincare.png");

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

ϑgrid = asm_grid(start=2, N=4)#, f_quad=1)
φgrid = asm_grid(start=-2, N=3)#, f_quad=1)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);

#bounds chosen to not go outside the bounding surfaces
sgrid = MID.ContGridDataT(N=50, start=0.4, stop=0.7)
cont_grids = init_grids(sgrid, ϑgrid, φgrid)

#%%
#now we construct the continuum
cont_ω = reconstruct_ω(MID.Spectrum.qfm_continuum(prob, cont_grids, surfs), cont_grids);
cont_norm_ω = reconstruct_ω(continuum(un_prob, cont_grids, false), cont_grids);


#%%
continuum_plot(cont_ω, cont_grids, filename="data/qfm/qfm_cont.png")
continuum_plot(cont_norm_ω, cont_grids, filename="data/qfm/unpert_cont.png")



#%%

#now we consider the full spectrum

#interesting, the left hand side seems to be cooked no matter what!
#running this from 0.4 (where the flux surfaces are good!) causes r=0.4-0.5 to be completly cooked
#but running this from 0.2 causes it to be cooked from 0.2 to 0.4, but perfect from 0.4-0.5. 
#Perhaps we need to change the island amplitude!
#with the amplitude quadratic, and a smaller amp, going from 0.4 seems to have fixed the issue, 
#Seems like the problemo is just that the flux surface at the left edge is both the largest, and least dense part of the flux surfaces
rgrid = MID.Structures.rfem_grid(N=80, start=0.4, stop=0.8)
θgrid = MID.Structures.asm_grid(start=2, N=4)#, f_quad=1)
ζgrid = MID.Structures.asm_grid(start=-2, N=3)#, f_quad=1)#, incr=2)

grids = init_grids(rgrid, θgrid, ζgrid)

#%%
#seems to be a fair bit slower!
#makes sense, will want to use MIDParallel oneday.
evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, surfs=surfs, full_spectrum=true);
evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, surfs=surfs, full_spectrum=false, nev=200, target_freq=0.3);

evals_norm, ϕ_norm, ϕft_norm = compute_spectrum(prob=un_prob, grids=grids, full_spectrum=true);

evals_weird, ϕ_weird, ϕft_weird = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);

#%%


continuum_plot(evals)#), filename="data/qfm/qfm_spectrum.png")
continuum_plot(evals_norm)#, filename="data/qfm/norm_spectrum.png")
continuum_plot(evals_weird, filename="data/qfm/weird_spectrum.png")

#=

continuum_plot(evals, filename="data/qfm/qfm_full_spectrum.png")
continuum_plot(evals_norm, filename="data/qfm/norm_full_spectrum.png")
continuum_plot(evals_weird, filename="data/qfm/weird_full_spectrum.png")
=#


#%%

ind_norm = find_ind(evals_norm, 0.1703)
ind = find_ind(evals, 0.1703)

potential_plot(ϕft_norm, grids, ind_norm, filename="data/qfm/norm_tae.png")
potential_plot(ϕft, grids, ind, filename="data/qfm/tae.png")
