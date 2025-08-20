#doing interpolation for the derivative
#ideally, we want to be able to compute the second derivative at a point.
#unsure if we need to interpolate a larger region
using MID
#%%

f(x) = x^2 * exp(x)
df(x) = 2*x*exp(x) + x^2 * exp(x)
ddf(x) = 2*exp(x) + 4*x*exp(x) + x^2 * exp(x)

N1 = 15
xg = LinRange(0, 1, N1)

fvals = zeros(N1, 2)

for i in 1:N1
    fvals[i, 1] = f(xg[i])
    fvals[i, 2] = df(xg[i])
end
#%%

function deriv_int(x, fvals, xgrid, order)

    ind = argmin(abs.(x .- xgrid))

    ξ, xinds, dx = MID.Mapping.global_to_local(ind, xgrid, x)

    grid_id = [0, 0, 1, 1]
    basis_id = [0, 1, 0, 1]
    val = 0.0
    for hx in 1:4
        gi = xinds[grid_id[hx] + 1]
        bi = 1 + basis_id[hx]
        val += fvals[gi, bi] * MID.Mapping.hb(ξ, hx, dx, order)
    end
    return val
end
deriv_int(0.7, fvals, xg, 1)
df(0.7)

df(1.0)
deriv_int(xg[2], fvals, xg, 2)
ddf(xg[2])
