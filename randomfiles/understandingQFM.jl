
using MID
using MIDViz
using Plots; plotlyjs()
using Plots; gr()
#try to understand the action function, from after the coefficeints are returned.
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
contourf(r2D_α
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
