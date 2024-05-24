
using MID
using Plots; plotlyjs() #having this here, and installed in the global environment tricks it into using plotlyjs for interactive plots. This is an awful solution.
#also gives some fkn warning, I think becuase MID doesn't have PlotlyJS.


#this file shows basis usage for common aspects of MID
#should also be a good verification to make sure everything is working as it should.
#doesn't have anything to do with islands or damping though!

N = 200;
#the collect is a bit annoying, but ok because we will typically use a clustered grid.
rgrid = collect(LinRange(0, 1, N));

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo); 
grids = init_grids(rgrid=rgrid, mstart=2, mcount=2, nstart=-2, ncount=1);
#tae_freq = (0.395845)^2/10^2

#inputs_to_file(prob=prob, grids=grids, dir="data/")

ω_cont = continuum(prob=prob, grids=grids);


plot_continuum(ω = ω_cont, rgrid=rgrid)


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true);

"og"
3×3 Matrix{Float64}:
 0.0  -0.0187637     0.000244251
 0.0  -0.0857724     0.000453065
 0.0   0.000490373  -2.49237e-6

 "New guy"
3×3 Matrix{Float64}:
  0.0          -0.00938184    0.000122125
 -0.00938184   -2.24294e-18   0.000471719
  0.000122125   0.000471719  -3.97047e-23


reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)

tae_ind = find_ind(ω, 0.380)
plot_potential(r=rgrid, ϕ=ϕ, ind=tae_ind, pmd=grids.pmd, n=1)

tae_freq = (ω[tae_ind] / geo.R0)^2

#now solve using Arpack with higher resolution.
N = 2000;
#the collect is a bit annoying, but ok because we will typically use a clustered grid.
rgrid = collect(LinRange(0, 1, N));

#problem stays the same.
grids = init_grids(rgrid=rgrid, mstart=2, mcount=2, nstart=-2, ncount=1);

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq);

reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)

plot_potential(r=rgrid, ϕ=ϕ, ind=1, pmd=grids.pmd, n=1)

#0.3803

a = [1, 2, 3]
b = [3, 5, 6]

using LinearAlgebra
dot(a, b)

