"""
Continuum Damping

In this example we compute the continuum damping of a TAE.
This example follows the problem outlined in Bowden and Hole 2015.
"""

using ChaoticShearAlfvenWaves
using CSAWViz

#this example uses an aspect ratio of 10.
geo = init_geometry(:tor, R0=10.0)

#The radial variable is set to the geometric radius with :r.
#we also have the specific q and density profile needed to 'close' the frequnecy gap 
fields = init_fields(:r, q=damping_q, dens=damping_dens)

#set the artificial damping parameter, δ, within the finite-Larmor-radius struct.
flr = init_flr(δ=-4.0e-9)

#combine to create the problem
prob = init_problem(geometry=geo, fields=fields, flr=flr)

# The TAE interacts with the continuum at roughly r≈0.765
#so the grid is clustered, with about 50% of the points around this region.
rgrid = init_grid(:r, 1000, sep1=0.75, sep2=0.78, frac=0.5)
θgrid = init_grid(:sm, 2, start=1)
ζgrid = init_grid(:sm, 1, start=-1)

grids = init_grids(rgrid, θgrid, ζgrid)

#target the TAE at 0.323, defaults to 100 eigenvalues.
solver = init_solver(prob=prob, target=0.323)

evals, ϕ, ϕft = compute_spectrum(prob, grids, solver);

#looking closer at the TAE
ind = find_ind(evals, 0.323)

#We know have a complex frequency, showing damping.
display(evals.ω[ind])
#the damping ratio
#note that this is different that predicted by Bowden and Hole
#as we solve the full SAW equation, rather than the reduced form.
display(imag(evals.ω[ind]) / real(evals.ω[ind])) 

#Harmonics show the strong interaction with the continuum.
harmonic_plot(ϕft, grids, ind)



