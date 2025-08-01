
using MID
using MIDCantori
using MIDViz
using Plots; plotlyjs()
using Plots; gr()
#try to understand the action function, from after the coefficeints are returned.
#%%

k = 0.003
geo = init_geo(R0=1.0)
isl1 = init_island(m0=4, n0=-3, A=k/4, flux=true)
isl2 = init_island(m0=3, n0=-2, A=k/3, flux=true)
isls = MID.IslandT[isl1, isl2]

prob = init_problem(q=cantori_q, isls=isls, geo=geo, type=:flux, met=:cylinder)

#%%
r0, θ0 = MIDCantori.QFM.anal_periodic_orbit(a, b, 2, k)

poincare_plot(prob, 1000, collect(LinRange(0.58, 0.62, 50)), zeros(50), xlimits=(0, 2π), ylimits=(0.59, 0.61))
poincare_plot(prob, 1000, [0.6], [0.0], xlimits=(0, 2π), ylimits=(0.59, 0.61))
scatter!(LinRange(0, 2π/a, length(νarr)), gl[1] .+ νarr)#, ylimits=(0.55, 0.65))

#why on earth is ν a constant for each psuedo field line, that seems to be actually wild
#wild that a fixed correction can create a periodic orbit.
#also seems odd that this is a perfect sin function,
#perhaps this can be analytically computed based on the size of perturbation?
#also looks like this could place upper limits on the qfm process.
gl = surface_guess([(13, 9)], cantori_q)
#%%

αvals = LinRange(0, 2π/a, length(νarr))
k*(sin(3*αvals[3]) + sin(4*αvals[3]))
k*(cos(3*αvals[3])/3 + cos(4*αvals[3])/4)

νarr[3] #* gl[1]^2

scatter(αvals, νarr)
scatter!(αvals, @. k*(cos(3*αvals)/3 + cos(4*αvals)/4))
#%%
a = 13
b = 9
r1 = (a, b)

met = MID.MetT()
B = MID.BFieldT()
M = 32
N = 8
rcosarr, rsinarr, θcosarr, θsinarr, νarr = MID.QFM.action(r1, prob, met, B, M, N, 0.7, 2);

scatter(LinRange(0, 2π/a, length(νarr)), νarr)
#%%
#size if (number of field lines, number of coefficients (a*N+1)).
size(rcosarr)

#first step is to convert back to normal space,
r = MID.QFM.irfft1D(rcosarr, rsinarr, 1)
θ = MID.QFM.irfft1D(θcosarr, θsinarr, 1)
ζ = LinRange(0, 2π*a, 2*N*a)

θf = zeros(8, 2*N*a)
for i in 1:8
    θf[i, :] = @. θ[i, :] + b/a * ζ
end


fll = 3
scatter(θ[fll, 1:2*N:end], r[fll, 1:2*N:end])
scatter(θf[fll, 1:2*N:end], r[fll, 1:2*N:end])

2π / (b) * 2
6π/5

2π * (b/a)
2π/b

mod(b*3, a)


LinRange(0, 2π/7, 9)[1:end-1]

θf[:, 1:2*N:end]
mod.(θf[:, 1:2*N:end], 2π)
2π/a
4π/a

#%%
p = plot()
for i in 1:8
    scatter!(θ[i, 1:2*N:end], r[i, 1:2*N:end])
end
display(p)


#so we have the solution in the domain (0, 2π/a), we need to extend the field lines to 2π.
#%%
r2dal = zeros((a*8, a*8));
θ2dal = zeros((a*8, a*8));
for i in 0:a-1
    idx = mod(b*i, a)
    display(idx)
    r2dal[1+idx*2*N:(idx+1)*2*N, 1:(a-i) * N*2] = r[:, 1+i*N*2:end]
    r2dal[1+idx*2*N:(idx+1)*2*N, 1+(a-i) * N*2: end] = r[:, 1:i*N*2]
    θ2dal[1+idx*2*N:(idx+1)*2*N, 1:(a-i) * N*2] = θ[:, 1+i*N*2:end]
    θ2dal[1+idx*2*N:(idx+1)*2*N, 1+(a-i) * N*2: end] = θ[:, 1:i*N*2]
end

#%%
fll = 4
#ok, so we only consider α [0, 2π/a] as the rest are the same
#this looks to duplicated each field line, but start from the next position essentially.
#but for some reason this is transposed....
scatter(θ2dal[fll, 1:2*N:end], r2dal[fll, 1:2*N:end])
scatter!(θ2dal[fll+8, 1:2*N:end], r2dal[fll+8, 1:2*N:end])
scatter(θ2dal[fll, 1:2*N:6*N], r2dal[fll, 1:2*N:6*N])
scatter!(θ2dal[fll+8, 1:2*N:6*N], r2dal[fll+8, 1:2*N:6*N])
scatter(θ2dal[fll, 1:N:end], r2dal[fll, 1:N:end])
scatter!(θ2dal[fll+8, 1:N:end], r2dal[fll+8, 1:N:end])
scatter(θ2dal[1:2*N:end, fll], r2dal[1:2*N:end, fll])
scatter(θ2dal[1:end, fll], r2dal[1:end, fll])

scatter(θ2dal[fll, 2*N:2*N], r2dal[fll, 2*N:2*N])
scatter!(θ2dal[fll, 4*N:4*N], r2dal[fll, 4*N:4*N])
scatter!(θ2dal[fll+8, 1:2*N:end], r2dal[fll+8, 1:2*N:end])
#%%
myr2d = zeros((8*a, 2*N*a));
myθ2d = zeros((8*a, 2*N*a));

#lets see if we can do this another way.
#%%
n_fl = 8
#doesn't work as the starting point needs to be shifted for each batch of the field lines.
#so the first fl always hits the same a points, but the order is different.
for i in 1:7
    display((i-1)*n_fl)
    display((i)*n_fl)
    myr2d[1+(i-1)*n_fl:(i*n_fl), :] = r[:, :]
    myθ2d[1+(i-1)*n_fl:(i*n_fl), :] = θ[:, :]
end

display(myr2d .- r2dal)

scatter(myθ2d[fll, 1:2*N:end], myr2d[fll, 1:2*N:end])
scatter!(myθ2d[fll+8, 1:2*N:end], myr2d[fll+8, 1:2*N:end])
scatter(myθ2d[fll, 1:2*N:6*N], myr2d[fll, 1:2*N:6*N])
scatter!(myθ2d[fll+8, 1:2*N:6*N], myr2d[fll+8, 1:2*N:6*N])
#%%
#ok now we have r(α, ζ), θ(α, ζ)
#need to change the α to ϑ via ϑ = α + b/a ζ

#%%

k1 = 0.005
geo = init_geo(R0=4.0)
isl1 = init_island(m0=3, n0=-2, A=k1/3)

isls = [isl1]#, isl2, isl3]

prob = init_problem(geo=geo, q=low_shear_qfm_q, isls=isls)
#%%
Ntraj = 150
#rlist = collect(LinRange(0.7, 0.95, Ntraj));
rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

poincare_plot(prob, Nlaps, Ntraj, rlist)
#%%
rats1 = lowest_rationals(7, prob.q(0.0)[1], prob.q(1.0)[1])
gl1 = surface_guess(rats1, prob.q)
@time surfs1 = construct_surfaces(rats1, gl1, prob, M=32, N=16);
plot_surfs(surfs1)

#%%
rats = [(3, 2), (4, 3)]
rats = [(5, 3)]
gl = surface_guess(rats, prob.q)
#rcosarr, rsinarr, θcosarr, θsinarr, νarr = construct_surfaces(rats, gl, prob, M=16, N=8)
met = MID.MetT()
B = MID.BFieldT()
rcosarr, rsinarr, θcosarr, θsinarr, νarr = MID.QFM.action(rats[1], prob, met, B, 16, 8, gl[1], 2);
#%%
αpoints, fp = size(rcosarr)
size(rsinarr)

ζarr = LinRange(0, 2*rats[1][1]*π, fp + 2*(floor(Int64, (fp-1)-fp/2)+1))[1:end-1]
r = zeros((αpoints, length(ζarr)));
θ = zeros((αpoints, length(ζarr)));
a, b = rats[1]
#b, a = rats[1]
for αind in 1:αpoints

    r[αind, :] .= rcosarr[αind, 1]
    θ[αind, :] .= θcosarr[αind, 1]
    
    for n in 1:fp-1
        r[αind, :] += @. rcosarr[αind, n+1] * cos(n*ζarr/a) + rsinarr[αind, n+1]*sin(n*ζarr/a)
        θ[αind, :] += @. θcosarr[αind, n+1] * cos(n*ζarr/a) + θsinarr[αind, n+1]*sin(n*ζarr/a)
    end
end
#%%
plot(ζarr, r[3, :])
plot(ζarr, θ[3, :])

#%%
#so there is some complex shit going to that prevent us doing this the same, probbaly doesn't matter
rft = MID.QFM.irfft1D(rcosarr, rsinarr);
θcosarr[:, 1] .= 0.0 #each theta curve starts at a different value based on the intial α value
#this condition seems to remove that!
ζarr = LinRange(0, 2*a*π, size(rft)[end]+1)[1:end-1];
θft = MID.QFM.irfft1D(θcosarr, θsinarr);
#%%
plot(ζarr, rft[3, :])
plot(ζarr, θft[3, :])
#%%
p = plot()
for i in 1:αpoints
#for i in 1:5
    plot!(θft[i, :])
end
display(p)

#%%
#lets see if we can figure out wot the fek r2d alpha is.
#the fact that these are the dimensions is already a fkn disaster
N = 8
MM = 4
fM = N * MM
r2D_α = zeros((αpoints*a, αpoints * a)) .+ 0.82;
θ2D_α = zeros((αpoints*a, αpoints * a));
#for i in 0:a-1
for i in 0:a-1
    idx = mod(b*i, a)
    r2D_α[1+idx*fM:(idx+1)*fM, 1:(a-i)*fM] = rft[:, 1+i*fM:end]
    #r2D_α[1+idx*fM:(idx+1)*fM, 1+(a-i)*fM:end] = rft[:, 1:i*fM]
    θ2D_α[1+idx*fM:(idx+1)*fM, 1:(a-i)*fM] = θft[:, 1+i*fM:end]
    #θ2D_α[1+idx*fM:(idx+1)*fM, 1+(a-i)*fM:end] = θft[:, 1:i*fM]
    #θ2D_α[1+idx*fM:(idx+1)*fM, 1+(a-i)*fM:end] = θft[:, 1+i*fM:end]
    #r2D_α[1+idx*fM:(idx+1)*fM, 1:(a-i)*fM] .= 0.0
end
#%%
contourf(r2D_α)
contourf(θ2D_α)

#%%
contourf(r[:, 1+1*32:end])
contourf(r[:, 1:1*32])
contourf(θft)
contourf(rft)
#%%
r2D_ϑ = zeros((αpoints*a, fM));
θ2D_ϑ = zeros((αpoints*a, fM));

for i in 0:fM-1
    arr = 0:a*fM
    idx = @.mod(arr-i*b, a*fM) .+ 1
    idx = idx[1:end-1]
    r2D_ϑ[:, i+1] = r2D_α[idx, i+1]
end

contourf(r2D_ϑ)
