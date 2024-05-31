"""

This file gives a basic use case of using the MID package to find TAE. This file does not consider a magnetic island chain or continuum damping. Please see island_damping.jl in examples for this.

This file first defines the grids and the problem, this is then solved for the full spectrum using Julia's inbuilt eigenvalue solver. This allows us to find the frequency of the TAE.

Then the grid resolution is increased and we solve using Arpack, using shift-and-invert to solve for only a few eigenvalues near the previously found TAE frequency.
"""


using MID
using Plots; plotlyjs() #having this here, and installed in the global environment tricks it into using plotlyjs for interactive plots. This is an awful solution.
#also gives some fkn warning, I think becuase MID doesn't have PlotlyJS.

#may want to use a more extreme example of a tae, eg fu_dam, so TAE is clearer.

#this file shows basis usage for common aspects of MID
#should also be a good verification to make sure everything is working as it should.
#doesn't have anything to do with islands or damping though!
"""

We start with the grids, we define the number of radial points to use.
Then we create the grids using init_grids(). In this simple case, a radial grid will be made from 0 to 1 with N points, we consider 2 poloidal modes, 2 and 3, while using a single toroidal mode n=-2.

"""

N = 100;
grids = init_grids(N=N, mstart=2, mcount=2, nstart=-2, ncount=1);

"""

Then we define the problem to be solved, first we define our geometry, which just contains the major radius. Then we define the problem, which requires a q-profile and the geometry. The problem will default to a uniform density and a toroidal metric.

""";

geo = GeoParamsT(R0=10.0)
prob = init_problem(q=Axel_q, geo=geo); 

"""

We can also write and read this setup from file if desired, which is useful when comparing cases to have all the inputs stored.

""";

#inputs_to_file(grids=grids, prob=prob, dir="data/examples/");

#prob, grids = inputs_from_file(dir="data/examples/");


"""

Now we can compute the continuum, this extract only the second radial derivative terms that contain the singularity defining the continuum. This does not work once island chains are introduced.

Compared to the full solver, this is much faster as each radial point is considered one at a time, reducing the size of the matrix to solve.

Followed by plotting the continuum, with the different modes labelled

""";

ω_cont = continuum(prob=prob, grids=grids);
plot_continuum(ω = ω_cont, grids=grids)

"""

Now we solve the full spectrum case, this returns the normalised eigenvalues and the eigenfunctions.
The flag full_spetrum, means the inbuild Julia solver is used, which is much slower than 
Arpack, but returns all of the eigenvalues.

"""

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true);


"""

From the eigenfunctions, we can crudely reconstruct the continuum, this is done by taking the location of the maximum value of the potential, and mapping it with the corresponding eigenvalue.

This is a flawed method, but retains most of the continuum structure and highlights any global modes inside gaps.

"""

reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)


"""

Next we find the index of the TAE, determining the frequency is most easily done with an interactive plotting library, such as PlotlyJs(). In this case the TAE occurs at ~0.380.

"""

tae_ind = find_ind(ω, 0.380)

"""

Now we can plot the mode structure of the tae

"""

plot_potential(grids=grids, ϕ=ϕ, ind=tae_ind, n=1)


"""

Now that we are certain we have found a TAE, we determine the un-normalised TAE frequency for the shift-and-invert algorithm.

"""

tae_freq = (ω[tae_ind] / geo.R0)^2

"""

Now we can employ a higher resolution radial grid, which becomes very important for accurate damping calculations.

We redefine the grids used, while the problem remains the same.

"""

N = 2000;
grids = init_grids(N=N, mstart=2, mcount=2, nstart=-2, ncount=1);


"""

Now we solve, noting that full_spectrum is set to false (default behaviour), and we have passed in the target frequency.

"""

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq);


"""

Again, we can reconstruct the continuum, but now we only have a few eigenvalues.
Then we can plot the mode structure again, noting that now the TAE is located at index 1.

"""

reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)
plot_potential(grids=grids, ϕ=ϕ, ind=1, n=1)

