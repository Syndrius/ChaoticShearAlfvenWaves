


Nr = 50;
Nθ = 5

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo); 

rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=2)
ζgrid = init_sm_grid(start=-2, count=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)



ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true);


#ϕms = mode_structure(ϕ, grids);
#reconstruct_continuum(ω, ϕms, grids)#, filename="data/fu_dam_spectrum.png")

#unsure why it is flipped... hopefully a resolution problemo
#this is an extremly different tae frequency... real good lol
tae_ind = find_ind(ω, 0.396)
tae_freq = ω[tae_ind]

@test tae_ind == 98
@test tae_freq ≈ 0.3963 atol=0.001

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq);


tae_ind = find_ind(ω, 0.396)
tae_freq = ω[tae_ind]

@test tae_ind == 1
@test tae_freq ≈ 0.3963 atol=0.001