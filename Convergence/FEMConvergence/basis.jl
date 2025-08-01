
#something is wrong with the basis, here we solve the 1d bessel equation to try and fix the relative amplitude of the 
#derivative part of the basis.
#this should allow us to interpolate properly
#hopefully this fixes a few tthings throughout the code.
using MID
using Plots
#%%


N = 10
xgrid = LinRange(0, π, N)
#xgrid = LinRange(0, 1, N)
gp = 5
llam, lf = MID.Helmholtz.linear(N, gp)

clam, cf = MID.Helmholtz.cubic(N, gp);

nlam, nf = MID.Helmholtz.no_integration(N, gp);
clam
nlam
#%%

clam[5]
plot(xgrid, cf[1:2:end, 5])
plot!(xgrid, cf[2:2:end, 5])
nlam[5]
plot!(xgrid, nf[1:2:end, 5])
plot!(xgrid, nf[2:2:end, 5])

plot!(xgrid, 0.33*sin.(3 .*xgrid))
plot!(xgrid, cos.(3 .*xgrid))

ξ, wg = gausslegendre(3)
S = MID.Basis.hermite_basis(ξ)
#%%
dx = xgrid[2] - xgrid[1]
tot = 0.0
for i in 1:3
    tot += S.H[1, i] * S.H[1, i] * wg[i] * (dx / 2)
end
display(tot)

13 / 35 * dx
    
