
#the interpolation package is fkn terrible
#so we need to get the hermite stuff working
#here we test in 1d to make sure it works properly
using MID
using Plots; gr()

#start simple
#f(x) = sin(x)
#df(x) = cos(x)
f(x) = x*(1-x)*exp(x)
df(x) = -exp(x)*(x^2+x-1)

N1 = 10
#xvals = LinRange(0.0, 2π, N1+1)[1:end-1]
xvals = LinRange(0.0, 1, N1)
fvals = zeros(N1, 2)
fvals[:, 1] = f.(xvals)
fvals[:, 2] = df.(xvals)
#fvals = f.(xvals)
#dfvals = df.(xvals)
#%%

plot(xvals, f.(xvals))
#%%
function hb(ξ, h, Δx)
    #additional jacobian term is very awkward.

    #local to global jacobain is no 2/ Δx
    #so we need to invert the local to global tangent scaling.

    if h==1
        return MID.Basis.h00(ξ) #/ Δx
    elseif h==2
        return MID.Basis.h10(ξ) * Δx / 2
    elseif h==3
        return MID.Basis.h01(ξ) #/ Δx
    else
        return MID.Basis.h11(ξ) * Δx / 2
    end
end
function global_to_local(ind, grid, x)
    #guard in case point is exactly on grid
    if grid[ind]==x
        #extra edge cases only for derivatives.
        ξ = -1.0
        ind1 = ind
        #assumes 2π periodicity!
        if ind1 == length(grid)
            #display("seeing periodicity?")
            ind2 = 1
            Δx = 2π + (grid[ind2] - grid[ind1])  #global difference
        else
            ind2 = ind+1 #doesn't matter
            Δx = grid[ind2] - grid[ind] 
        end
        inds = [ind1, ind2]
    #periodic case
    elseif x > grid[end]
        inds = [length(grid), 1]
        Δx = 2π + (grid[1] - grid[end])
        #ξ = (x - grid[end]) / dx
        #this should be the only difference with new basis
        ξ = (x - grid[end]) * 2 / Δx - 1
    elseif grid[ind] < x
        inds = [ind, ind+1]
        Δx = grid[ind+1] - grid[ind]
        #ξ = (x - grid[ind]) / Δx
        ξ = (x - grid[ind]) * 2 / Δx - 1
    else
        #may need to add some other edge cases 
        #in particular, the qfm grid not being maximal can cause problemos
        #typically this is easily fixed by making the mapped grid smaller.
        inds = [ind-1, ind]
        Δx = grid[ind] - grid[ind-1]
        ξ = (x - grid[ind-1]) * 2 / Δx - 1
        #ξ = (x - grid[ind-1]) / Δx
    end
    return ξ, inds, Δx

end
#%%
function int(x, fvals, xvals)
    ind = argmin(abs.(x .- xvals))
    #display(abs.(x .- xvals))

    ξ, inds, Δx = global_to_local(ind, xvals, x)
    #grid_id = [0, 0, 1, 1]
    #basis_id = [0, 1, 0, 1]
    grid_id = MID.Indexing.grid_id
    basis_id = MID.Indexing.basis_id

    val = 0.0
    for hx in 1:4
        gi = inds[grid_id[hx]+1]
        bi = 1 + basis_id[hx]
        val += fvals[gi, bi] * hb(ξ, hx, Δx)
    end

    return val
end
display(int(0.67, fvals, xvals))
display(f(0.67))
#%%
N2 = 100
#xi = LinRange(0.0, 2π, N2+1)[1:end-1]
xi = LinRange(0.0, 1, N2)
ϕ = zeros(N2);

for (i, x) in enumerate(xi)
    #ϕ[i] = int(x, fvals, xvals)
    ϕ[i] = MID.Mapping.hermite_interpolation(x, fvals .+ 0.0im, collect(xvals))
end

#%%

plot(xi, ϕ)
plot!(xvals, fvals[:, 1])
plot!(xi, f.(xi))
#plot!(xvals, fvals)

#%%
#now lets consider the 2d case.
f(x, y) = x^2*sin(y) #only periodic in 1d now.
dfdx(x, y) = 2*x*sin(y)
dfdy(x, y) = x^2*cos(y)
d2fdxdy(x, y) = 2*x*cos(y)

N1 = 30
xvals = LinRange(0.0, 1, N1)
yvals = LinRange(0.0, 2π, N1+1)[1:end-1]
fvals = zeros(N1, N1, 4);
fvals[:, :, 1] = f.(xvals, yvals')
fvals[:, :, 2] = dfdx.(xvals, yvals')
fvals[:, :, 3] = dfdy.(xvals, yvals')
fvals[:, :, 4] = d2fdxdy.(xvals, yvals')

#dfdxvals = dfdx.(xvals, yvals')
#dfdyvals = dfdy.(xvals, yvals')
#d2fdxdyvals = d2fdxdy.(xvals, yvals')
#%%

function int_2d(x, y, fvals, xvals, yvals)
    xind = argmin(abs.(x .- xvals))
    yind = argmin(abs.(y .- yvals))

    grid_id = [0, 0, 1, 1]
    basis_id = [0, 1, 0, 1]

    ξx, xinds, Δx = global_to_local(xind, xvals, x)
    ξy, yinds, Δy = global_to_local(yind, yvals, y)
    #grid_id = [0, 0, 1, 1]
    #basis_id = [0, 1, 0, 1]
    grid_id = MID.Indexing.grid_id
    basis_id = MID.Indexing.basis_id

    val = 0.0
    for hx in 1:4, hy in 1:4
        gi = (xinds[grid_id[hx]+1], yinds[grid_id[hy]+1])

        #most likely needs to be swapped as our structure here is different to ϕ in the normal case.
        #yep, this is in the opposit eorder to ϕ, no big deal!
        bi = 1 + basis_id[hx] + 2*basis_id[hy]
        val += fvals[gi..., bi] * hb(ξx, hx, Δx) * hb(ξy, hy, Δy)
    end
    return val
end
#%%
display(int_2d(0.3, 0.5, fvals, xvals, yvals))
#println(collect(xvals))
#println(collect(yvals))
display(f(0.3, 0.5))
#%%

N2 = 200
xi = LinRange(0.0, 1.0, N2)
yi = LinRange(0.0, 2π, N2+1)[1:end-1]
ϕ = zeros(N2, N2);

for (i, x) in enumerate(xi), (j, y) in enumerate(yi)
    ϕ[i, j] = int_2d(x, y, fvals, xvals, yvals)
end
#%%

contourf(ϕ)
contourf(f.(xi, yi'))



fi = f.(xi, yi')
contourf(ϕ .- fi)
maximum(abs.(ϕ .- fi))
#%%
plot(xi, fi[:, 25])
plot!(xi, ϕ[:, 25])
#%%
#now the finally, 3d!
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
function int_3d(x, y, z, fvals, xvals, yvals, zvals)
    xind = argmin(abs.(x .- xvals))
    yind = argmin(abs.(y .- yvals))
    zind = argmin(abs.(z .- zvals))
    #exact point!
    #just need to worry about the edges now, start by considering a grid that is only 
    #between the edges.
    if xvals[xind] == x
        ξx = 0.0
        xind1 = xind
        xind2 = xind #doesn't matter
        dx = 0.0 #shouldn't matter I think
    #x is not periodic now.
    elseif xvals[xind] < x

        #will need to check for endpoints!
        xind1 = xind
        xind2 = xind + 1

        dx = xvals[xind2] - xvals[xind1] #global difference
        ξx = (x - xvals[xind1]) / dx
        #ξ = (xvals[ind2] - x) / dx
    else
        xind1 = xind-1
        xind2 = xind
        dx = xvals[xind2] - xvals[xind1] #global difference
        ξx = (x - xvals[xind1]) / dx
        #ξ = (xvals[ind2] - x) / dx
    end
    if yvals[yind] == y
        ξy = 0.0
        yind1 = yind
        yind2 = yind #doesn't matter
        dy = 0.0 #shouldn't matter I think

    #periodic case
    elseif y > yvals[end]
        yind1 = length(yvals)
        yind2 = 1
        dy = 2π + (yvals[yind2] - yvals[yind1])  #global difference
        ξy = (y - yvals[yind1]) / dy

    elseif yvals[yind] < y

        #will need to check for endpoints!
        yind1 = yind
        yind2 = yind + 1

        dy = yvals[yind2] - yvals[yind1] #global difference
        ξy = (y - yvals[yind1]) / dy
        #ξ = (xvals[ind2] - x) / dx
    else
        yind1 = yind-1
        yind2 = yind
        dy = yvals[yind2] - yvals[yind1] #global difference
        ξy = (y - yvals[yind1]) / dy
        #ξ = (xvals[ind2] - x) / dx
    end
    if zvals[zind] == z
        ξz = 0.0
        zind1 = zind
        zind2 = zind #doesn't matter
        dz = 0.0 #shouldn't matter I think

    #periodic case
    elseif z > zvals[end]
        zind1 = length(zvals)
        zind2 = 1
        dz = 2π + (zvals[zind2] - zvals[zind1])  #global difference
        ξz = (z - zvals[zind1]) / dz

    elseif zvals[zind] < z

        #will need to check for endpoints!
        zind1 = zind
        zind2 = zind + 1

        dz = zvals[zind2] - zvals[zind1] #global difference
        ξz = (z - zvals[zind1]) / dz
        #ξ = (xvals[ind2] - x) / dx
    else
        zind1 = zind-1
        zind2 = zind
        dz = zvals[zind2] - zvals[zind1] #global difference
        ξz = (z - zvals[zind1]) / dz
        #ξ = (xvals[ind2] - x) / dx
    end

    grid_id = [0, 0, 1, 1]
    basis_id = [0, 1, 0, 1]

    xinds = [xind1, xind2]
    yinds = [yind1, yind2]
    zinds = [zind1, zind2]

    val = 0.0
    for hx in 1:4, hy in 1:4, hz in 1:4
        gi = (xinds[grid_id[hx]+1], yinds[grid_id[hy]+1], zinds[grid_id[hz]+1])

        #this lil trick doesn't work with our setup, hopefully it would actually work in the general case
        bi = 1 + 1*basis_id[hx] + 2*basis_id[hy] + 4*basis_id[hz]
        val += fvals[gi..., bi] * hb(ξx, hx, dx) * hb(ξy, hy, dy) * hb(ξz, hz, dz)
    end

    return val
end
#%%

display(int_3d(0.3, 0.5, 0.1, fvals, xvals, yvals, zvals))

display(f(0.3, 0.5, 0.1))

#%%
Ni1 = 200
Ni2 = 100
Ni3 = 20
xi = LinRange(0.0, 1.0, Ni1)
yi = LinRange(0.0, 2π, Ni2+1)[1:end-1]
zi = LinRange(0.0, 2π, Ni3+1)[1:end-1]
ϕ = zeros(Ni1, Ni2, Ni3);
fi = zeros(Ni1, Ni2, Ni3);

#why is this so quick?
for (i, x) in enumerate(xi), (j, y) in enumerate(yi), (k, z) in enumerate(zi)
    ϕ[i, j, k] = int_3d(x, y, z, fvals, xvals, yvals, zvals)
    fi[i, j, k] = f(x, y, z)
end
#%%
#coolio hoolio.
#just need to be v careful of the ol loops
contourf(ϕ[100, :, :])
contourf(fi[100, :, :])
#%%
#speed comparison

#this can only do linear, guess that is why it is quicker
f_itp = extrapolate(interpolate((xvals, yvals, zvals), fvals[:, :, :, 1], (Gridded(Cubic()), Gridded(Cubic()), Gridded(Cubic()))), Smooth() );


for (i, x) in enumerate(xi), (j, y) in enumerate(yi), (k, z) in enumerate(zi)
    ϕ[i, j, k] = int_3d(x, y, z, fvals, xvals, yvals, zvals)
end

#ok so noticibly quicker.
for (i, x) in enumerate(xi), (j, y) in enumerate(yi), (k, z) in enumerate(zi)
    ϕ[i, j, k] = f_itp(x, y, z)
end
