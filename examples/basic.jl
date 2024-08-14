"""
This is now basic case for 1d fem!! needs to be adjusted appropriatly.

This file gives a basic use case of using the MID package to find TAE. This file does not consider a magnetic island chain or continuum damping. Please see island_damping.jl in examples for this.

This file first defines the grids and the problem, this is then solved for the full spectrum using Julia's inbuilt eigenvalue solver. This allows us to find the frequency of the TAE.

Then the grid resolution is increased and we solve using Arpack, using shift-and-invert to solve for only a few eigenvalues near the previously found TAE frequency.
""";


using MID
using Plots; plotlyjs() #having this here, and installed in the global environment tricks it into using plotlyjs for interactive plots. This is an awful solution.
#also gives some fkn warning, I think becuase MID doesn't have PlotlyJS.
"""

We start with the grids. Each 1d grid can be defined for the finite element method or the spectral method. Here we consider finite elements in r and spectral in θ, ζ. This is the fastest method, useful for simple serial cases and to demonstrate how this package works.

Each grid is defined individually, the finite element grid requires the number of points:

""";

Nr = 100;
rgrid = init_fem_grid(N=Nr);

"""

Then we define the spectral grids, these require a starting mode then the total count of modes. Here we consider two poloidal modes and a single toroidal mode.

""";

θgrid = init_sm_grid(start=1, count=2);
ζgrid = init_sm_grid(start=-1, count=1);

"""

Now we combine the three grids into a single grid structure, which is one of the main inputs for our key functions.

""";

grids = init_grids(rgrid, θgrid, ζgrid);

"""

Then we define the problem to be solved, first we define our geometry, which just contains the major radius. For the Fu Van Damm case, we have R0=4.0, with a minor radius a=1 the default behaviour.

""";

geo = GeoParamsT(R0=4.0);

"""

Then we define the problem, which, at a minimum, requires a q-profile and the geometry. Further inputs are shown in the other example files. The problem will default to a uniform density and a toroidal metric. 

Some q-profiles are predefined, including q= 1 + r^2 used for this case. Defining other q-profiles is straightforward.

""";

q_profile = fu_dam_q;
prob = init_problem(q=q_profile, geo=geo); 

"""

With the grids and the problem structs defined, we now have all required inputs and can compute the spectrum. 

""";

#should chaneg cr to evals. Cause it is just evals but with some extra data.
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);


"""

This returns:
- evals, a struct containing eigenvalues and associated information on each eiegenvalue such as radial location and dominant fourier mode.
- ϕ, the eigenfunctions, which is a 4d matrix, the first index denotes which eigenfunction, then the remaining 3 are the 3d grid (r, θ, ζ).
- ϕft, the eigenfunctions but with a fourier transform in θ and ζ to easily look at the mode structure. This is in the same form as ϕ.

We can now view the reconstruction of the continuum:

""";

plot_continuum(evals);

"""

Now we can investigate a specific eigenvalue in more depth. First we find the index corresponding the eigenvalue, in this case the TAE.

""";

tae_ind = find_ind(evals, 0.29)
tae_freq = evals.ω[tae_ind]

"""

Now we plot the mode structure of the eigenfunction.

""";

plot_potential(ϕft, grids, tae_ind);

""";

Now that we know the frequency of interest, we can increase the grid resolution, target the specific region of the spectrum using a shift and invert algorithm to improve efficiency.

""";

Nr = 1000;
rgrid = init_fem_grid(N=Nr);
θgrid = init_sm_grid(start=1, count=2)
ζgrid = init_sm_grid(start=-1, count=1);
grids = init_grids(rgrid, θgrid, ζgrid)
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq);


#now conitnuum only shows region around the tae frequency
plot_continuum(evals);
#Tae is now the first index returned.
plot_potential(ϕft, grids, 1);

"""

Alternatively, we can employ finite elements in θ, noting that this method is more demanding computationally.

We also define a 'phase factor' for the θ grid, to focus on the mode with m=1

""";

Nr = 100;
Nθ = 10;
rgrid = init_fem_grid(N=Nr);
θgrid = init_fem_grid(N=Nθ, pf=1)
ζgrid = init_sm_grid(start=-1, count=1);
grids = init_grids(rgrid, θgrid, ζgrid)
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq);

plot_continuum(evals);

#in this low resolution case the tae frequency is slightly higher, but the methods converge with larger resolutions.
tae_ind = find_ind(evals, 0.30)
tae_freq = evals.ω[tae_ind]
plot_potential(ϕft, grids, tae_ind)

"""

Or we can consider finite elements in every dimension, however, this is slow to construct the matrices and requires parallel implementation to be of practical use.

""";

Nr = 50;
Nθ = 5;
Nζ = 2;
rgrid = init_fem_grid(N=Nr);
θgrid = init_fem_grid(N=Nθ, pf=1)
ζgrid = init_fem_grid(N=Nζ, pf=-1)
grids = init_grids(rgrid, θgrid, ζgrid)
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq);

plot_continuum(evals);

#in this low resolution case the tae frequency is slightly higher, but the methods converge with larger resolutions.
tae_ind = find_ind(evals, 0.326)
tae_freq = evals.ω[tae_ind]
#we are still able to find a tae, but the mode structure has not adequately resolved. 
plot_potential(ϕft, grids, tae_ind)


"""

Finally, we can also save our input structures and results to files, useful when larger cases are considered.

""";


Nr = 1000;
rgrid = init_fem_grid(N=Nr);
θgrid = init_sm_grid(start=1, count=2)
ζgrid = init_sm_grid(start=-1, count=1);
grids = init_grids(rgrid, θgrid, ζgrid)

geo = GeoParamsT(R0=4.0);
q_profile = fu_dam_q;
prob = init_problem(q=q_profile, geo=geo); 

#input structures are written to file
inputs_to_file(prob=prob, grids=grids, dir="data/examples/")

#then can then be read and run as normal.
prob_file, grids_file = inputs_from_file(dir="data/examples/");

#Alternatively, we can directly compute the spectrum from the file
spectrum_from_file(dir="data/examples/", σ=real.(tae_freq))

#this saves the output to the same directory.
evals = evals_from_file(dir = "data/examples/");

plot_continuum(evals);

#And we can read eigenfunctions from file. When run from file each eigenfunction is written individually.

ϕft = efunc_from_file(dir = "data/examples/", ind=1)
plot_potential(ϕft, grids)

