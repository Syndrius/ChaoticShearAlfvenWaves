#think at the end of this one, we might mention that other FLR things exist
#but are still experimental.

#Current set up adequatly replicated Bowden Singular.
#the frequency is a bit off due to difference in equations.
using MID
using Plots; plotlyjs()
using MIDViz
#%%

#how the fek are we initialising this????
#we want to specify) the geometry and the fields to match coordinates right?
#think the met should be given as a symbol, :torus, :cyl, :island etc.
geo = init_geometry(:tor, R0=10.0)

fields = init_fields(:r, q=damping_q, dens=damping_dens)

#does this need to be negative?
flr = init_flr(δ=-4.0e-9)
#%%
#perhaps within init_prob we actually set the met function properly?
#I am assuming we will have to recreate geo and fields
#which is probably ok.
#maybe we can just make a few placeholder metric and island
#functions.
#think that would work,
#so geo has the met and met_func, once of which will be set to whatever is passed in
#and the other is set as a placeholder. which is then updated later.
#may be a bit annoying that prob.geo won't match geo
#could cause problems if we want to pass in geo to stuff?
#all the same with fields tbh.
#feel like this is a very niche problem that we will probably ignore.
#but does it mean that the met and fields functions should actually be defined within prob
#and geo and fields initialisation is just for iitialisation?
#maybe that is a neater way of doing it?
#fields and geo contain the info, but not in a useful form
#and prob collates and instantiates it?
#think that would be better tbh!
#or is is better to have a Fields with different types 
#and just multi dipatch on that, passing fields directly into compute_B!()?
#seems like one way works well for geo and the other way works well for fields...
#might be good to have symbol storing the type for each of the individual structs
#then in problem we can assert that they are the same, or try and convert them
#if it is not possible we can error out.
#note that having an different geo etc is always going to be needed as the island metric is instantiated with an island, which is within fields...
#so geo will probably need to store 2 symbols:
#one for the actual geometry, :cyl, :tor etc
#one for flux vs radius vs kappa. #may even be able to incorporate qfm in here somehow.
prob = init_problem(geometry=geo, fields=fields, flr=flr)
#%%

rgrid = init_grid(:r, 1000, sep1=0.75, sep2=0.78, frac=0.5)
θgrid = init_grid(:sm, 2, start=1)
ζgrid = init_grid(:sm, 1, start=-1)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%
solver = init_solver(prob=prob, full_spectrum=true)
solver = init_solver(prob=prob, target=0.323)
#%%
#looks to be working now, just converging to a slightly different value
#seem to need 1000 points before the tae is clear!
evals, ϕ, ϕft = compute_spectrum(prob, grids, solver);

continuum_plot(evals)

#find_ind is no fkn good.
ind = find_ind(evals, 0.323)
ind = 1
evals.ω[ind]
imag(evals.ω[ind]) / real(evals.ω[ind]) 

#not perfect, but this example is adequate I think.
harmonic_plot(ϕft, grids, ind)



