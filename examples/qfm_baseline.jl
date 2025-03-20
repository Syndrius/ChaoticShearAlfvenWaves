
#example showing qfm for simple perturbation that changes the flux surfaces from striaght to slightly wibbly.

#%%
#so we need this now, perhaps worth making it default/always use by putting it in the startup
using Revise
using MID

using Plots; plotlyjs()
using JLD2

function qfm_bench_q(r::Float64)

    a = 1.0
    b = 2.0

    return a + b * r^2, 2 * b * r
end



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
#for comparisons!
prob2 = init_problem(q=qfm_bench_q, geo=geo, met=cylindrical_metric!);
#%%


#find the edge surfaces for constructing the Farey tree.
#bound at (1, 2), r~0.7
#and at (3, 4), r~0.4
@time bounding_surfs = construct_surfaces([1, 3], [2, 4], [0.4, 0.7], prob);

#looks like good bounding surfaces!
plot_surfs(bounding_surfs);
#%%


qlist, plist = farey_tree(5, 1, 2, 3, 4)

println(qlist)
println(plist)

#o it is not about these values, it is about when the q-profile equals these values... -> no, the iota profile equals these values.
#may want to introduce some kind of check?
println(qlist ./ plist)

guess_list = 0.5 .* ones(length(qlist))
#this seems to be cooked now, may need to compare with vscode?
#Ran for like ~25 mins without finishing?
@time surfs = construct_surfaces(qlist, plist, guess_list, prob);

#not really an even distribution of surfaces, but oh well.
plot_surfs(surfs);
#%%


#not certain if we want to save these and the itp or only one or the other?
save_object("data/benchmark_surfs.jld2", surfs)

#%%

surfs = load_object("data/benchmark_surfs.jld2");
#plot_surfs(surfs);

#%%

#surf_itp = MID.QFM.create_surf_itp(surfs);

#hmm so fquad may need to be higher than we expect.
ϑgrid = asm_grid(start=2, N=4)#, f_quad=1)
φgrid = asm_grid(start=-2, N=3)#, f_quad=1)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);

#bounds chosen to not go outside the bounding surfaces
sgrid = MID.ContGridDataT(N=50, start=0.42, stop=0.7)
grids = init_grids(sgrid, ϑgrid, φgrid)

#%%

#surf interpolation is created inside hte function
#given creating the interpolation is instant, probably fine.
#still all masive, but the transform function is cooked af.
ω = MID.Spectrum.qfm_continuum(prob, grids, surfs)

#%%

function reconstruct_ω(ω, grids)
    #turn a weir doutput thing into a normal output hopefully.
    new_ω = zeros(grids.r.N, grids.θ.N, grids.ζ.N)

    for i in 1:grids.r.N, j in 1:grids.θ.N * grids.ζ.N

        θ, ζ = MID.Indexing.index_to_grid(j, grids)
        new_ω[i, θ, ζ] = ω[i, j]
    end
    return new_ω

end

new_ω = reconstruct_ω(ω, grids)
continuum_plot(new_ω , grids, ymax=5)
#%%

#shape of ω is not right for this. I guess we need to reconstruct.
#continuum_plot(ω, grids)
rgrid = rfem_grid(N=50, start=0.42, stop=0.7)
θgrid = asm_grid(start=2, N=2)#, f_quad=1)
ζgrid = asm_grid(start=-2, N=1)#, f_quad=1)#, incr=2)

grids = init_grids(rgrid, θgrid, ζgrid)


#seems to be a fair bit slower!
#makes sense, will want to use MIDParallel oneday.
evals, ϕ, ϕft = compute_spectrum_qfm(prob=prob, grids=grids, surfs=surfs, full_spectrum=true);

#%%
continuum_plot(evals, ymax=1)

potential_plot(ϕft, grids, 80)


#%%

evals_norm, ϕ_norm, ϕft_norm = compute_spectrum(prob=prob2, grids=grids, full_spectrum=true);

continuum_plot(evals_norm)

potential_plot(ϕft_norm, grids, 10)

#%%
evals_weird, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);

continuum_plot(evals_weird)
#%%


rgrid = rfem_grid(N=50, start=0.42, stop=0.7)
θgrid = afem_grid(N=8)#, f_quad=5)
ζgrid = asm_grid(start=-2, N=2)#, f_quad=5)#, incr=2)

grids_ffs = init_grids(rgrid, θgrid, ζgrid)

#so this gives almost exactly the same garbage. in some ways that is a good thing!
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids_ffs, surfs=surfs, full_spectrum=true);

continuum_plot(evals, ymax=10^4)


#interest, clearly something is going wrong
#but this kind of looks like we are not ages off.
#look, promising start, but something is very wrong with this
#hopefully the continuum case is the same problemo as the other case.s
#%%
potential_plot(ϕ, grids_ffs, 300)
#%%

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
#Comparison to unpertrubed continuum



prob2 = init_problem(q=qfm_bench_q, geo=geo, met=cylindrical_metric!)

#hmm so fquad may need to be higher than we expect.
ϑgrid = asm_grid(start=-2, N=4)#, f_quad=1)
φgrid = asm_grid(start=-2, N=4)#, f_quad=1)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);

#bounds chosen to not go outside the bounding surfaces
sgrid = MID.ContGridDataT(N=50, start=0.42, stop=0.7)
grids = init_grids(sgrid, ϑgrid, φgrid)



#surf interpolation is created inside hte function
#given creating the interpolation is instant, probably fine.
#still all masive, but the transform function is cooked af.
ω_norm = reconstruct_ω(continuum(prob2, grids, false), grids)

ω_qfm = reconstruct_ω(MID.Spectrum.qfm_continuum(prob, grids, surfs), grids)

#%%

#seems to be pretty much identical!
##interesting
#Perhaps means our perturbation is a bit small, unfort this will require creating more surfaces...
#Which is cooked atm for some reason!
continuum_plot(ω_qfm , grids, ymax=5)
continuum_plot(ω_norm , grids, ymax=5)
#%%



#given creating the interpolation is instant, probably fine.
#still all masive, but the transform function is cooked af.

#%%

function reconstruct_ω(ω, grids)
    #turn a weir doutput thing into a normal output hopefully.
    new_ω = zeros(grids.r.N, grids.θ.N, grids.ζ.N)

    for i in 1:grids.r.N, j in 1:grids.θ.N * grids.ζ.N

        θ, ζ = MID.Indexing.index_to_grid(j, grids)
        new_ω[i, θ, ζ] = ω[i, j]
    end
    return new_ω

end

new_ω = reconstruct_ω(ω, grids)
continuum_plot(new_ω , grids, ymax=5)
#%%
