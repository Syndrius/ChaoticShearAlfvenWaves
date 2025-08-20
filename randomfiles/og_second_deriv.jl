
#using our interpolations to compute the second derivative
#for determining mode sharpness.
using Plots
using Plots; plotlyjs()

#start simple
f(x) = sin(x)
df(x) = cos(x)
ddf(x) = -sin(x)

N1 = 100
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



function dh00(t)

    return 6*t^2 - 6*t
    return 2*t^3 - 3*t^2 + 1
    #return (1+2*t)*(1-t)^2
end

function dh10(t)
    return (1-t)^2 - 2*t*(1-t)
    #return 2 * (t^3-2*t^2+t)
    return t*(1-t)^2
end

function dh01(t)
    return -6*t^2 + 6*t
    return -2t^3+3t^2
    #return t^2*(3-2*t)
end

function dh11(t)
    return 2*t*(t-1) + t^2
    #return 2*(t^3-t^2)
    return t^2*(t-1)
end

function ddh00(t)

    return 12*t - 6
    return 6*t^2 - 6*t
end

function ddh10(t)
    return -4*(1-t)+2*t
    return (1-t)^2 - 2*t*(1-t)
end

function ddh01(t)
    return -12*t + 6
    return -6*t^2 + 6*t
end

function ddh11(t)
    return 4*t + 2*(t-1)
    return 2*t*(t-1) + t^2
end

function hb(t::Float64, h::Int64, dt::Float64, deriv::Int64)
    #additional jacobian term is very awkward.

    if h==1
        if deriv == 1
            return dh00(t) / dt
        elseif deriv == 2
            return ddh00(t) / dt^2
        else
            return h00(t)
        end
    elseif h==2
        if deriv == 1
            return dh10(t) 
        elseif deriv == 2
            return ddh10(t) / dt
        else
            return h10(t) * dt
        end
    elseif h==3
        if deriv == 1
            return dh01(t) / dt
        elseif deriv == 2
            return ddh01(t) / dt^2
        else
            return h01(t)
        end
    else
        if deriv == 1
            return dh11(t) 
        elseif deriv == 2
            return ddh11(t) / dt
        else
            return h11(t) * dt
        end
    end
end
#%%

tvals = collect(0:0.01:1)
plot(tvals, h00.(tvals))
plot!(tvals, h10.(tvals))
plot!(tvals, h01.(tvals))
plot!(tvals, h11.(tvals))
#%%
function global_to_local(ind, grid, x)
    
    if grid[ind] == x
        ξ = 0.0
        ind1 = ind
        #assumes 2π periodicity!
        if ind1 == length(grid)
            ind2 = 1
            dx = 2π + (grid[ind2] - grid[ind1])  #global difference
        else
            ind2 = ind+1 #doesn't matter
            dx = grid[ind2] - grid[ind] 
        end
    #periodic case
    elseif x > grid[end]
        ind1 = length(grid)
        ind2 = 1
        dx = 2π + (grid[ind2] - grid[ind1])  #global difference
        ξ = (x - grid[ind1]) / dx
    elseif grid[ind] < x
        #will need to check for endpoints!
        ind1 = ind
        ind2 = ind + 1
        dx = grid[ind2] - grid[ind1] #global difference
        ξ = (x - grid[ind1]) / dx
    else
        ind1 = ind-1
        ind2 = ind
        dx = grid[ind2] - grid[ind1] #global difference
        ξ = (x - grid[ind1]) / dx
    end
    return ξ, [ind1, ind2], dx
end
#%%
function int(x, fvals, xvals, order)
    ind = argmin(abs.(x .- xvals))

    ξ, xinds, dx = global_to_local(ind, xvals, x)

    grid_id = [0, 0, 1, 1]
    basis_id = [0, 1, 0, 1]
    val = 0.0
    for hx in 1:4
        gi = xinds[grid_id[hx]+1]
        bi = 1 + basis_id[hx]
        val += fvals[gi, bi] * hb(ξ, hx, dx, order)
    end
    return val
end
println(collect(xvals))
display(int(2π-0.1, fvals, xvals, 2))
display(ddf(2π-0.1))
#%%
N2 = 200
xi = LinRange(0.0, 2π, N2+1)[1:end-1]
ϕ = zeros(N2);

for (i, x) in enumerate(xi)
    ϕ[i] = int(x, fvals, xvals, 2)
end

#%%

plot(xi, ϕ)
plot!(xi, ddf.(xi))
#plot!(xvals, fvals)
#%%
#now higher dim.
f(x, y, z) = x^2*sin(y)*cos(z) 
dfdx(x, y, z) = 2*x*sin(y)*cos(z)
dfdy(x, y, z) = x^2*cos(y)*cos(z)
dfdz(x, y, z) = -x^2*sin(y)*sin(z)
dfdxdy(x, y, z) = 2*x*cos(y)*cos(z)
dfdxdz(x, y, z) = -2*x*sin(y)*sin(z)
dfdydz(x, y, z) = -x^2*cos(y)*sin(z)
dfdxdydz(x, y, z) = -2*x*cos(y)*sin(z)
dfdx2(x, y, z) = 2*sin(y) * cos(z)
dfdx2dy(x, y, z) = 2*cos(y) * cos(z)
dfdx2dydz(x, y, z) = -2*cos(y) * sin(z)
dfdx2dz(x, y, z) = -2*sin(y) * sin(z)
dfdy2(x, y, z) = -x^2*sin(y)*cos(z)
dfdxdy2(x, y, z) = -2*x*sin(y)*cos(z)
dfdxdy2dz(x, y, z) = 2*x*sin(y)*sin(z)

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
function int_3d(x, y, z, fvals, xvals, yvals, zvals, xd, yd, zd)
    xind = argmin(abs.(x .- xvals))
    yind = argmin(abs.(y .- yvals))
    zind = argmin(abs.(z .- zvals))

    ξx, xinds, dx = global_to_local(xind, xvals, x)  
    ξy, yinds, dy = global_to_local(yind, yvals, y)  
    ξz, zinds, dz = global_to_local(zind, zvals, z)  
    grid_id = [0, 0, 1, 1]
    basis_id = [0, 1, 0, 1]

    val = 0.0
    for hx in 1:4, hy in 1:4, hz in 1:4
        gi = (xinds[grid_id[hx]+1], yinds[grid_id[hy]+1], zinds[grid_id[hz]+1])

        #this lil trick doesn't work with our setup, hopefully it would actually work in the general case
        bi = 1 + 1*basis_id[hx] + 2*basis_id[hy] + 4*basis_id[hz]
        val += fvals[gi..., bi] * hb(ξx, hx, dx, xd) * hb(ξy, hy, dy, yd) * hb(ξz, hz, dz, zd)
    end

    return val
end
#%%

display(int_3d(0.3, 0.5, 0.1, fvals, xvals, yvals, zvals, 2, 1, 1))

display(dfdx2dydz(0.3, 0.5, 0.1))

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
    ϕ[i, j, k] = int_3d(x, y, z, fvals, xvals, yvals, zvals, 2, 1, 1)
    fi[i, j, k] = dfdx2dydz(x, y, z)
end
#%%
#coolio hoolio.
#looks to be working pretty well.
contourf(ϕ[100, :, :])
contourf(fi[100, :, :])
contourf(ϕ[100, :, :] .- fi[100, :, :])
#%%
