#note that this one has been updated for newest compute_B!
#other cases have not!
#maybe the convergence testing should be in test?


using MID
using Plots; plotlyjs()
#this file is for B/vectors rather than matrices
#this function will stand in for the different aspects we want to compute
#J and magB look good
function f(r::Float64, θ::Float64, ζ::Float64, prob::MID.ProblemT)

    met = MID.MetT() #predefining causes issues with multiple calls at the same time.

    MID.Geometry.island_metric!(met, r, θ, ζ, prob.geo.R0, prob.isls[1])
    return met.gu
end

function df(r, θ, ζ, prob)

    met = MID.MetT()
    #MID.Equilibrium.compute_B!(B, met, prob.q, prob.isls, r, θ, ζ)
    MID.Geometry.island_metric!(met, r, θ, ζ, prob.geo.R0, prob.isls[1])

    #return B.dmag_B
    return met.dgu
end
function dfdr(r, θ, ζ, prob, Δr)
    #display(f(r+Δr, θ, ζ, prob))
    return (f(r+Δr, θ, ζ, prob) - f(r - Δr, θ, ζ, prob)) / (2 * Δr)
    return (-f(r+2Δr, θ, ζ, prob) + 8 * f(r + Δr, θ, ζ, prob) - 8 * f(r-Δr, θ, ζ, prob) + f(r - 2*Δr, θ, ζ, prob)) / (12 * Δr)
end
function dfdθ(r, θ, ζ, prob, Δθ)
    #display("where and why")
    return (f(r, θ+Δθ, ζ, prob) - f(r, θ-Δθ, ζ, prob)) / (2 * Δθ)
end
function dfdζ(r, θ, ζ, prob, Δζ)
    return (f(r, θ, ζ+Δζ, prob) - f(r, θ, ζ-Δζ, prob)) / (2 * Δζ)
end
#%%
geo = init_geo(R0=1.0)
isl1 = init_island(m0=4, n0=-3, A=0.005, coords=true, qp=2.0, ψ0=0.5)
prob = init_problem(geo=geo, q=island_q, isl=isl1, type=:island)
#%%
hlist = [0.01, 0.005, 0.003, 0.001, 0.0005, 0.0003, 0.0001]
Nh = length(hlist)
Nr = 13
Nθ = 13
Nζ = 13
rvals = LinRange(0.2, 0.8, Nr)
θvals = LinRange(0.2, 6.0, Nθ)
ζvals = LinRange(0.2, 6.0, Nζ)

error = zeros(Nh, Nr, Nθ, Nζ, 3, 3, 3);

for (i, h) in enumerate(hlist), (j, r) in enumerate(rvals), (k, θ) in enumerate(θvals), (l, ζ) in enumerate(ζvals)
    an = df(r, θ, ζ, prob)
    #display((r, θ, ζ))

    error[i, j, k, l, 1, :, :] = abs.(an[:, :, 1] .- dfdr(r, θ, ζ, prob, h))
    error[i, j, k, l, 2, :, :] = abs.(an[:, :, 2] .- dfdθ(r, θ, ζ, prob, h))
    error[i, j, k, l, 3, :, :] = abs.(an[:, :, 3] .- dfdζ(r, θ, ζ, prob, h))
end

avg_error = zeros(Nh, 3, 3, 3);

for i in 1:Nh, j in 1:3, k in 1:3, l in 1:3
    avg_error[i, j, k, l] = maximum(error[i, :, :, :, j, k, l])
end
lh = log.(hlist)
le = log.(avg_error)

p = plot()#ylims=(0, 0.1))
for i in 1:3, j in 1:3, k in 1:3
    #plot!(-lh, le[:, i, j, k], label="d"*string(i)*"g_"*string(j) * string(k))
    plot!(-lh, avg_error[:, i, j, k], label="d"*string(i)*"g_"*string(j) * string(k))
end
display(p)


#%%
grad = zeros(Nh-1, 3, 3, 3)

for i in 1:Nh-1
    grad[i, :, :, :] = @. (le[i+1, :, :, :] - le[i, :, :, :]) / (lh[i+1] - lh[i])
end
p = plot()
for i in 1:3, j in 1:3, k in 1:3
    #plot!(-lh, le[:, i, j, k])
    plot!(-lh[1:end-1], grad[:, i, j, k], label="d"*string(i)*"g_"*string(j) * string(k))
end
display(p)



#%%
#yikes, something is wrong with dgu...

function tor_test(r, θ, ζ)

    #peak efficiency
    met = MID.Geometry.MetT()
    MID.Geometry.toroidal_metric!(met, r, θ, ζ, 4.0)
    t1 = inv(met.gl)

    #display(met.gu .- t1)

    test = zeros(3, 3, 3)
    for i in 1:3
        test[:, :, i] = -t1 * met.dgl[:, :, i] * t1
    end
    #display(met.dgu .- test)

    display(met.gu)
    display(test[:, :, 1])
    display(met.dgu[:, :, 2])
    display(met.dgu[:, :, 3])
end

tor_test(0.1, 0.0, 0.0)


#%%

dgl = df(0.4, 0.4, 0.4)[:, :, 1]
fdgl = dfdr(0.4, 0.4, 0.4, 0.3)
log.(abs.(dgl .- fdgl))
