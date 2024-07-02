

#start very small, matrix scales much more extremly
Nr = 100;

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo)#, met=no_delta_metric!); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=2, count=2)
ζgrid = init_sm_grid(start=-2, count=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)



ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true);


#reconstruct_continuum(ω, ϕ, grids)#, filename="data/fu_dam_spectrum.png")

#unsure why it is flipped... hopefully a resolution problemo
tae_ind = find_ind(ω, 0.39)
tae_freq = ω[tae_ind]

@test tae_ind == 198
@test tae_freq ≈ 0.3764 atol=0.001


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ = tae_freq);


tae_ind = find_ind(ω, 0.39)
tae_freq = ω[tae_ind]

@test tae_ind == 1
@test tae_freq ≈ 0.3764 atol=0.001