
#the interpolation package is fkn terrible
#so we need to get the hermite stuff working
#here we test in 1d to make sure it works properly
using Plots

#start simple
f(x) = sin(x)
df(x) = cos(x)

N1 = 20
xvals = LinRange(0.0, 2π, N1+1)[1:end-1]
fvals = zeros(N1, 2)
fvals[:, 1] = f.(xvals)
fvals[:, 2] = df.(xvals)
#fvals = f.(xvals)
#dfvals = df.(xvals)
#%%

plot(xvals, f.(xvals))
#%%

function h00(t)

    return 2*t^3 - 3*t^2 + 1
    #return (1+2*t)*(1-t)^2
end

function h10(t)
    #return 2 * (t^3-2*t^2+t)
    return t*(1-t)^2
end

function h01(t)
    return -2t^3+3t^2
    #return t^2*(3-2*t)
end

function h11(t)
    #return 2*(t^3-t^2)
    return t^2*(t-1)
end

function hb(t, h, dt)
    #additional jacobian term is very awkward.

    if h==1
        return h00(t)
    elseif h==2
        return h10(t) * dt
    elseif h==3
        return h01(t)
    else
        return h11(t) * dt
    end
end
#%%

tvals = collect(0:0.01:1)
plot(tvals, h00.(tvals))
plot!(tvals, h10.(tvals))
plot!(tvals, h01.(tvals))
plot!(tvals, h11.(tvals))
#%%
function int(x, fvals, dfvals, xvals)
    ind = argmin(abs.(x .- xvals))
    #display(abs.(x .- xvals))

    #exact point!
    #just need to worry about the edges now, start by considering a grid that is only 
    #between the edges.
    if xvals[ind] == x
        #yet to test this case occuring!
        ξ = 0.0
        ind1 = ind
        ind2 = ind #doesn't matter
        #occuring!
        #display("occuring")
        dx = 0.0 #shouldn't matter I think

    #periodic case
    elseif x > xvals[end]
        ind1 = length(xvals)
        ind2 = 1
        dx = 2π + (xvals[ind2] - xvals[ind1])  #global difference
        ξ = (x - xvals[ind1]) / dx

    elseif xvals[ind] < x

        #will need to check for endpoints!
        ind1 = ind
        ind2 = ind + 1

        dx = xvals[ind2] - xvals[ind1] #global difference
        ξ = (x - xvals[ind1]) / dx
        #ξ = (xvals[ind2] - x) / dx
    else
        ind1 = ind-1
        ind2 = ind
        dx = xvals[ind2] - xvals[ind1] #global difference
        ξ = (x - xvals[ind1]) / dx
        #ξ = (xvals[ind2] - x) / dx
    end
    #display(xvals[ind1])

    #guesss, but gets the gist across
    #we will also need to generalise this to a loop over the values, perhaps one for ind1 and one for ind2? we will need to create an S object I think.

    grid_id = [0, 0, 1, 1]
    basis_id = [0, 1, 0, 1]

    xinds = [ind1, ind2]
    val = 0.0
    for hx in 1:4

        gi = xinds[grid_id[hx]+1]
        bi = 1 + basis_id[hx]
        val += fvals[gi, bi] * hb(ξ, hx, dx)
    end

    #display(x)
    #display(ξ)
    #ind1 denotes the ind of the grid point to the left of the x value
    #ind2 denotes the ind of the grid point to the right of the x value
    #need a dx as we are changing coords to local, so we get a jacobian term
    #makes this perf
    #val = fvals[ind1] * h00(ξ) + dfvals[ind1] * h10(ξ) * dx + fvals[ind2]*h01(ξ) + dfvals[ind2]*h11(ξ) * dx
    #without the deriv it seems to work pretty well
    #but the deriv is cooked af.
    #val = fvals[ind2] * h01(ξ) + fvals[ind1] * h00(ξ)
    #this approxes the deriv very well
    #are we not supposed to use the extra shape functions?
    #surely we are??
    #val = dfvals[ind2] * h01(ξ) + dfvals[ind1] * h00(ξ)
    #val = dfvals[ind2] * h11(ξ) + dfvals[ind1] * h10(ξ)

    return val
end
println(collect(xvals))
display(int(2π-0.1, fvals, dfvals, xvals))
display(f(2π-0.1))
#%%
N2 = 200
xi = LinRange(0.0, 2π, N2+1)[1:end-1]
ϕ = zeros(N2);

for (i, x) in enumerate(xi)
    ϕ[i] = int(x, fvals, dfvals, xvals)
end

#%%

plot(xi, ϕ)
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
    #we will want some clarity about the fkn division stuff.
    #it is a global_to_local jacobian transformation.
    #naturally we will want these functions within Basis.
    #val = fvals[ind1] * h00(ξ) + dfvals[ind1] * h10(ξ) * dx + fvals[ind2]*h01(ξ) + dfvals[ind2]*h11(ξ) * dx

    #display(xind)
    #display(ξx)
    #display(yind)
    #display(ξy)

    grid_id = [0, 0, 1, 1]
    basis_id = [0, 1, 0, 1]

    xinds = [xind1, xind2]
    yinds = [yind1, yind2]

    val = 0.0
    for hx in 1:4, hy in 1:4
        gi = (xinds[grid_id[hx]+1], yinds[grid_id[hy]+1])

        #most likely needs to be swapped as our structure here is different to ϕ in the normal case.
        #yep, this is in the opposit eorder to ϕ, no big deal!
        bi = 1 + basis_id[hx] + 2*basis_id[hy]
        val += fvals[gi..., bi] * hb(ξx, hx, dx) * hb(ξy, hy, dy)
    end

    #this was the difference, need to consider all 4 corners of the square.
    #gonna get messy real quick in 3d.
    #need to loopify this.
    #val = fvals[xind1, yind1] * h00(ξx)*h00(ξy) + fvals[xind2, yind2] * h01(ξx)*h01(ξy) + fvals[xind1, yind2] * h00(ξx)*h01(ξy) + fvals[xind2, yind1] * h01(ξx) * h00(ξy)

    #val += (dfdxvals[xind1, yind1] * h10(ξx)*h00(ξy) + dfdxvals[xind2, yind2] * h11(ξx)*h01(ξy) + dfdxvals[xind1, yind2] * h10(ξx)*h01(ξy) + dfdxvals[xind2, yind1] * h11(ξx) * h00(ξy)) * dx

    #val += (dfdyvals[xind1, yind1] * h00(ξx)*h10(ξy) + dfdyvals[xind2, yind2] * h01(ξx)*h11(ξy) + dfdyvals[xind1, yind2] * h00(ξx)*h11(ξy) + dfdyvals[xind2, yind1] * h01(ξx) * h10(ξy)) * dy

    #val += (d2fdxdyvals[xind1, yind1] * h10(ξx)*h10(ξy) + d2fdxdyvals[xind2, yind2] * h11(ξx)*h11(ξy) + d2fdxdyvals[xind1, yind2] * h10(ξx)*h11(ξy) + d2fdxdyvals[xind2, yind1] * h11(ξx) * h10(ξy)) * (dx * dy)

    #val += (dfdxvals[xind1, yind1] * h10(ξx)*h00(ξy) + dfdxvals[xind2, yind2] * h11(ξx)*h01(ξy)) *dx

    #val += (dfdyvals[xind1, yind1] * h00(ξx)*h10(ξy) + dfdyvals[xind2, yind2] * h01(ξx)*h11(ξy)) * dy

    #val += (d2fdxdyvals[xind1, yind1] * h10(ξx)*h10(ξy) + d2fdxdyvals[xind2, yind2] * h11(ξx)*h11(ξy)) * (dx * dy)

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
