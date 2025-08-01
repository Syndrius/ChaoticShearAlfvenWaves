using MID
using Plots
#%%
Nx = 8
Ny = 3
Nz = 3


#Cartesian case is at least working now, 
#so we can finally run this thorugh MID
#this will allow us to check fff vs fss etc,
#and hopefully check the phase factors
#still would be good to get the cylindrical version going
#so that we can be more confident for a less straightforward case
#however, it is still probably worth running this on Gadi to see how it looks.
clam, cf = MID.Helmholtz.cubic(Nx, Ny, Nz, 3, 3, 3);
clam[46]
#%%
xgrid = LinRange(0, 1, Nx)
ygrid = LinRange(0, 2π, Ny+1)[1:end-1]
zgrid = LinRange(0, 2π, Nz+1)[1:end-1]

plot(xgrid, real.(cf[46, :, 3, 3]))
plot(ygrid, real.(cf[2, 5, :, 3]))
plot(zgrid, real.(cf[2, 5, 3, :]))

#hard to tell if this are ok,
#but the seem fine.
contourf(ygrid, xgrid, real.(cf[45, :, :, 2]))
contourf(zgrid, xgrid, real.(cf[45, :, 2, :]))
contourf(zgrid, ygrid, real.(cf[1, 5, :, :]))

#%%
bc1 = sort(MID.Helmholtz.compute_cubic_bcs(7, 5))
bc2 = sort(MID.Helmholtz.cubic_bcs_cart(7, 5))
sort(bc2) .- sort(bc1)
