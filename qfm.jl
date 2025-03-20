
#%%
using Revise
using MID
using MIDViz
using JLD2

using Plots; gr()#plotlyjs()

function qfm_benchmark_q(r::Float64)
    #a = 0.954545
    #b = 1.515151
    
    a = 1.93333
    b = 1.66666

    return a + b*r^2, 2 * b * r
end

#%%
#define the problem to solve

#should be arbitrary!
R0=10.0

#amp needs further thought!
#define the non-resonant island
isl = IslandT(m0=3, n0=2, A=0.012)

geo = GeoParamsT(R0=R0)

prob = init_problem(q=qfm_benchmark_q, met=cylindrical_metric!, geo=geo, isl=isl)
un_prob = init_problem(q=qfm_benchmark_q, met=cylindrical_metric!, geo=geo)


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
qlist, plist = farey_tree(4, 2, 1, 3, 1)


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
sgrid = MID.ContGridDataT(N=50, start=0.2, stop=0.8)
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

rgrid = rfem_grid(N=80, start=0.2, stop=0.8)
θgrid = asm_grid(start=2, N=4)#, f_quad=1)
ζgrid = asm_grid(start=-2, N=3)#, f_quad=1)#, incr=2)

grids = init_grids(rgrid, θgrid, ζgrid)

#%%
#seems to be a fair bit slower!
#makes sense, will want to use MIDParallel oneday.
evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, surfs=surfs, full_spectrum=true);

evals_norm, ϕ_norm, ϕft_norm = compute_spectrum(prob=un_prob, grids=grids, full_spectrum=true);

evals_weird, ϕ_weird, ϕft_weird = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);

#%%

continuum_plot(evals, filename="data/qfm/qfm_spectrum.png")
continuum_plot(evals_norm, filename="data/qfm/norm_spectrum.png")
continuum_plot(evals_weird, filename="data/qfm/weird_spectrum.png")

