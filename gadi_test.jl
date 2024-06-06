
#don't think we have access to github, so need to use scp, so modifications here won't carry over

#so repl now works on Gadi, and we have confirmed that MID will run as expected.
#TODO
#Make sure we can submit a job with MID, ie check installation location of julia etc! this works!!!
#Make sure MIDParallel will work in login node #this works, now using cpardiso as the main solver as that is intel and pre-installed!
#get MIDParallel to work in job.

using MID
#using Plots; plotlyjs() #having this here, and installed in the global environment tricks it into using plotlyjs for interactive plots. This is an awful solution.
#also gives some fkn warning, I think becuase MID doesn't have PlotlyJS.


#this file shows basis usage for common aspects of MID
#should also be a good verification to make sure everything is working as it should.
#doesn't have anything to do with islands or damping though!

N = 100;
#the collect is a bit annoying, but ok because we will typically use a clustered grid.
#rgrid = collect(LinRange(0, 1, N));

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo); 
grids = init_grids(N=N, mstart=2, mcount=2, nstart=-2, ncount=1);
#tae_freq = (0.395845)^2/10^2

#inputs_to_file(prob=prob, grids=grids, dir="data/")

#ω_cont = continuum(prob=prob, grids=grids);


#plot_continuum(ω = ω_cont, rgrid=rgrid)


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true);

display(ω[198])

#reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)

#tae_ind = find_ind(ω, 0.380)
#plot_potential(r=rgrid, ϕ=ϕ, ind=tae_ind, pmd=grids.pmd, n=1)
