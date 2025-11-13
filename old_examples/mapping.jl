
#testinig the island to 
using MID

κgrid = init_grid(type=:rf, N=100, stop=0.999)
ᾱgrid = init_grid(type=:af, N=30)
τgrid = init_grid(type=:af, N=15)
isl_grids = init_grids(κgrid, ᾱgrid, τgrid)
dir = "/Users/matt/phd/MIDParallel/data/mapping/"

#%%
#naturally, this is slow af.
#this is perhaps a bit to slow, I guess it is the interpolation that is the problemo?
#we may need to try and make an interpolation object storing the polynomials or something, 
#ok so interpolation is defs the problemo. Need to construct an interpolation object!
#our method is unnacetably slow. Unsure how to fix tbh.
#perhaps we use the normal interpolation for testing, then just use our shit version
#once we have good results? And just cop the huge time sink
#periodic stuff will be annoying af though
@profview MID.PostProcessing.tor_spectrum_to_isl(dir, isl_grids)
