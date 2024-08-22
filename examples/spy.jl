
#have a look at the sparse matrices

using MID
#using UnicodePlots
using Plots#; plotlyjs()
using SparseArrays


#4x5x3 probably shows the pattern the best, but we will want to plot this withot java probably for real discussion.
rgrid = init_fem_grid(N=5);
θgrid = init_fem_grid(N=5, pf=2);
#θgrid = init_sm_grid(start=1, count=4)
ζgrid = init_fem_grid(N=5, pf=-2);
#ζgrid = init_sm_grid(start=1, count=3)

grids = init_grids(rgrid, θgrid, ζgrid);

geo = GeoParamsT(R0=10.0)

isl = IslandT(m0=5, n0=4, A=0.0)
prob = init_problem(q=Axel_q, geo=geo, isl=isl); 


#with @views. 22.907080 seconds (8.07 M allocations: 721.379 MiB, 1.36% gc time)
#fk load more allocations and gc without views.
#outrageous spead up shifting the ϕ[:, test, :, ...] to ϕ[:, testr, testθ, :, :]


W, I = construct(prob, grids);

#will be Nr blocks, each block contins 32 rows (8 * Nθ * Nζ), and 96 cols (3 * 32), 3 comes from the overlap. (except the first and last block which have less columns, which will be ignored! and boundary blocks which will also be ignored.)
#regardless of the split, the diagonal will always be full, (unless diagonal is larger than the blocks, which is unlikely in general case)

#so each row will have 3 * 8 * Nθ * Nζ non-zeros. We just need to assume that the diagonal takes as many as it can, and any left over go in to the non-diagonal.
#128-32

#128+96

#160-(128-32)
#findnz(W);

#findall(!iszero, W)

#p = scatter(ylimits=(0, 8*25))
#spy!(real.(W))

#awful way to zoom but does seem to work for some reason...
#this will be useful for our write up on anpreallocation
#which we should do...
display(xlims!(ylims!(spy(real.(W), legend=false), (205,395)), (205,395)))

spy(real.(W), legend=false, ylims=(0, 8*25), dpi=600)

savefig("CanZoom.png")