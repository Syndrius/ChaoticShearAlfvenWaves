
using MID

#so I guess the hermite stuff is fine, just the cost of doing cubic interpolation instead of linear...
#we may need to come up with a way to create a qfm grid that will be islandy. we will probably just settle for super high density in the island region.
#%%
f(x, y, z) = x^2*sin(y)*cos(z) #only periodic in 1d now.
dfdx(x, y, z) = 2*x*sin(y)*cos(z)
dfdy(x, y, z) = x^2*cos(y)*cos(z)
dfdz(x, y, z) = -x^2*sin(y)*sin(z)
dfdxdy(x, y, z) = 2*x*cos(y)*cos(z)
dfdxdz(x, y, z) = -2*x*sin(y)*sin(z)
dfdydz(x, y, z) = -x^2*cos(y)*sin(z)
dfdxdydz(x, y, z) = -2*x*cos(y)*sin(z)

N1 = 100
N2 = 50
N3 = 10
xvals = LinRange(0.0, 1, N1)
yvals = LinRange(0.0, 2π, N2+1)[1:end-1]
zvals = LinRange(0.0, 2π, N3+1)[1:end-1]
fvals = zeros(N1, N2, N3, 8);
for (i, x) in enumerate(xvals), (j, y) in enumerate(yvals), (k, z) in enumerate(zvals)
    fvals[i, j, k, 1] = f(x, y, z)
    fvals[i, j, k, 2] = dfdx(x, y, z)
    fvals[i, j, k, 3] = dfdy(x, y, z)
    #these need to be swapped for the loop to work, not ideal! Surely can work in real case though
    #as these kind of loops are how we construct the damn thing.
    fvals[i, j, k, 5] = dfdz(x, y, z)
    fvals[i, j, k, 4] = dfdxdy(x, y, z)
    fvals[i, j, k, 6] = dfdxdz(x, y, z)
    fvals[i, j, k, 7] = dfdydz(x, y, z)
    fvals[i, j, k, 8] = dfdxdydz(x, y, z)
end
#%%
#fvals is our stand in for ϕ.
#first just checking if this is noticibly slower than the other case.
Ni1 = 200
Ni2 = 100
Ni3 = 20
xi = LinRange(0.0, 1.0, Ni1)
yi = LinRange(0.0, 2π, Ni2+1)[1:end-1]
zi = LinRange(0.0, 2π, Ni3+1)[1:end-1]
ϕ = zeros(ComplexF64, Ni1, Ni2, Ni3);
fvals = fvals .+ 0.1im;
for (i, x) in enumerate(xi), (j, y) in enumerate(yi), (k, z) in enumerate(zi)
    ϕ[i, j, k] = MID.Mapping.hermite_interpolation(x, y, z, fvals, collect(xvals), yvals, zvals)
end
#perhaps this wasn't actually that slow, just a bit slower than the actual interpolation
#so I guess we can use this then?
