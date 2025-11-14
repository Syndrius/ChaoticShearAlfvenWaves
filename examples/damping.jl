
#replicate Axel's damping case.
using MID
#%%

#how the fek are we initialising this????
#we want to specify) the geometry and the fields to match coordinates right?
#think the met should be given as a symbol, :torus, :cyl, :island etc.
geo = init_geometry(10.0, met=MID.Geometry.rad_toroidal_metric!)

fields = init_fields(q=axel_q, dens=axel_dens)

#does this need to be negative?
flr = init_flr(δ=1e-7)
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
prob = init_problem(geo=geo, fields=fields, flr=flr)
#%%

rgrid = init_grid(:r, 100, sep1=0.85, sep2=0.95, frac=0.5)
θgrid = init_grid(:sm, 2, start=2)
ζgrid = init_grid(:sm, 1, start=-2)

grids = init_grids(rgrid, θgrid, ζgrid)


