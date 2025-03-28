#here we are testing our implementation of metric derivatives.
#start with the basic case, then we will compare the qfm surfaces stuff.
#shown that the derivs for our normal B and met are fine.
#still need to test qfm stuff.
using FiniteDifferences
using MID
using JLD2

f(x) = x^2
df(x) = 2*x

xvals = 0.0:0.1:2

yvals = f.(xvals)

dyvals = df.(xvals)

adyvals  = ForwardDiff.derivative.(f, xvals)
#%%
surfs = load_object("data/qfm/benchmark_surfs.jld2");
plot_surfs(surfs)
surf_itp = MID.QFM.create_surf_itp(surfs);

#ok so now we can do some serious work
#%%
function tor_met(x)

    r, θ, ζ = x

    #don't care about efficient atm

    CT = MID.QFM.CoordTsfmT()
    tor_met = MID.Geometry.MetT()
    qfm_met = MID.Geometry.MetT()
    tor_B = MID.Equilibrium.BFieldT()
    qfm_B = MID.Equilibrium.BFieldT()

    R0 = 100.0
    isl1 = init_island(m0=3, n0=2, A=0.005)
    isl2 = init_island(m0=1, n0=-4, A=0.0000)

    MID.QFM.coord_transform!(r, θ, ζ, CT, surf_itp)
    #MID.Geometry.toroidal_metric!(tor_met, r, θ, ζ, R0)
    MID.Geometry.toroidal_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], R0)
    MID.Equilibrium.compute_B!(tor_B, tor_met, fu_dam_q, isl1, isl2, CT.coords[1], CT.coords[2], CT.coords[3])
    MID.QFM.met_transform!(tor_met, qfm_met, CT)

    MID.QFM.B_transform!(tor_B, qfm_B, qfm_met, CT)
    #MID.Equilibrium.compute_B!(B, met, fu_dam_q, isl1, isl2, r, θ, ζ)

    #dgl seems to be good!
    #dgu also seems to be good! v nice.
    #main problemo is probably more the qfm surfaces tbh
    #dB is also no problemo!
    #mag_b is also good! v nice.
    #CT seems to be good. Still need to actually check the new metric and B.
    #display(qfm_met.gl .- (transpose(CT.JM) * tor_met.gl * CT.JM))
    return qfm_B.B[3]
end


function dtor_met(r, θ, ζ)

    #don't care about efficient atm

    #met = MID.Geometry.MetT()
    CT = MID.QFM.CoordTsfmT()
    tor_met = MID.Geometry.MetT()
    qfm_met = MID.Geometry.MetT()
    tor_B = MID.Equilibrium.BFieldT()
    qfm_B = MID.Equilibrium.BFieldT()

    R0 = 100.0
    isl1 = init_island(m0=3, n0=2, A=0.005)
    isl2 = init_island(m0=1, n0=-4, A=0.0000)

    MID.QFM.coord_transform!(r, θ, ζ, CT, surf_itp)
    #MID.Geometry.toroidal_metric!(tor_met, r, θ, ζ, R0)
    MID.Geometry.toroidal_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], R0)
    MID.Equilibrium.compute_B!(tor_B, tor_met, fu_dam_q, isl1, isl2, CT.coords[1], CT.coords[2], CT.coords[3])
    MID.QFM.met_transform!(tor_met, qfm_met, CT)

    MID.QFM.B_transform!(tor_B, qfm_B, qfm_met, CT)

    #MID.Equilibrium.compute_B!(B, met, fu_dam_q, isl1, isl2, r, θ, ζ)
    #MID.Geometry.cylindrical_metric!(met, r, θ, ζ, R0)

    return qfm_B.dB[3, :]
end
#%%

#so looks like the qfm_dgl is a bit wrong, seems to be concentrated in the φ direction
#wonder if this is just because tor_dgl is relativly large in the ζ components, because they have R0 factors?
#also unsure how much this matters.

#damn daniel the CT stuff is extremely slow, wot the hek
rvals = LinRange(0.4, 0.7, 10)

θvals = LinRange(0.001, 2π, 10)
ζvals = LinRange(0.001, 2π, 10)

#test = [0.6, 4.887144, 0.69990205]
#(9, 3, 10) is a huge different for dB[3, 3]. This seems to be the most problomatic component in general.
#issue seems to go away with large R0 values, I think this is just becuase all the vars are extremely small.
#%%
test = [rvals[9], θvals[4], ζvals[10]]
fd = grad(central_fdm(5, 1), tor_met, test)[1]
tor = dtor_met(test...)
#%%

display((abs.(fd .- tor) ./ tor))


err_threshold = 0.05
tol = 1e-4
#note dgu[2, 2] fails for 1e-10, this is just because the vallues are v large and there is some rounding errors
#passes for 1e-8
for (i, r) in enumerate(rvals), (j, θ) in enumerate(θvals), (k, ζ) in enumerate(ζvals)
    fd = grad(central_fdm(5, 1), tor_met, [r, θ, ζ])[1] 
    tor = dtor_met(r, θ, ζ)
    diff = abs.(fd .- tor) 
    if diff[1] > tol
        diff[1] = diff[1] / tor[1]
    end
    if diff[3] > tol
        diff[3] = diff[3] / tor[3]
    end
    if diff[2] > tol
        diff[2] = diff[2] / tor[2]
    end
    #relative diff doesn't really work for very small values.
 
    if diff[1] > err_threshold
        display((i, j,  k))
    elseif diff[2] > err_threshold
        display((i, j,  k))
    elseif diff[3] > err_threshold
        display((i, j,  k))
    end
end

