
using MID
using MIDViz
using Plots; plotlyjs()
#try to understand the action function, from after the coefficeints are returned.
#%%

MM = 4
M = 24
N = 8
met = MID.Geometry.MetT()
B = MID.Equilibrium.BFieldT()
R0=4.0

#amp needs further thought!
#define the non-resonant island
k = 0.0035
isl = init_island(m0=5, n0=-2, A=k/5)
#isl2 = init_island(m0=7, n0=-3, A=k/7)

geo = init_geo(R0=R0)

#to solve non-Hermitian
#flr = MID.Structures.FLRT(δ = 1e-18)
prob = init_problem(q=qfm_benchmark_q, geo=geo, isl=isl)#, isl2=isl2)#, flr=flr)

qlist, plist = farey_tree(4, 2, 1, 3, 1)

#%%
#ok so , the rcosarr etc, are 2d arrays, [i, j], the i index gives a different solution of the gradient, in other words it defines a specific field line that satisfied the action. Each one differs on their starting value of $θ$, which is equivalent to the area (unclear).
#the j index tells us the values of the ρ_n, which for the non-resonant case is a bit dull, as mainly the n=0 index is going off. correspondnig to j=8
#with a resonant island, the j=2 case is going off. The indexing is not super clear, but probably not that important.
p = 2
q = 5
rcosarr, rsinarr, θcosarr, θsinarr, nvarr = MID.QFM.action2(p, q, prob, met, B, MM, M, N);
#%%
#not too complicated I think, turns the fourier spikes into a smooth function
#each i index still denotes a different trajectory, which are all out of phase.
r = MID.QFM.irfft1D(rcosarr, rsinarr);
ζ = LinRange(0, 2*q*π, size(r)[end]+1)[1:end-1];
θcosarr[:, 1] .= 0.0;
θ = MID.QFM.irfft1D(θcosarr, θsinarr);

#%%
fM = MM * N
qfM = q * fM
Nfft= qfM

r2D_alpha = zeros((qfM, Nfft));
θ2D_alpha = zeros((qfM, Nfft));

#this indexing is cooked af.
for i in 0:q-1
    idx = mod(p*i, q)

    #wot the actual fuck is going on here.
    r2D_alpha[1+idx * fM : (idx+1)*fM, 1 :(q-i) * N * MM] = r[:, 1+ i * N * MM : end]

    r2D_alpha[1+idx * fM : (idx+1)*fM, 1 + (q-i) * N * MM : end] = r[:, 1: i * N * MM]

    θ2D_alpha[1+idx * fM : (idx+1)*fM, 1 :(q-i) * N * MM] = θ[:, 1+i * N * MM : end]

    θ2D_alpha[1+idx * fM : (idx+1)*fM, 1 + (q-i) * N * MM : end] = θ[:, 1: i * N * MM]
end
#%%

#unclear what this is actually doing. Should be striaghtening the usrfaces??
#result of this step is not clear, I guess it probably shouldn't looks like anything? as this is r as a function of ϑ etc, so probably should looks rubbish

r2D_vartheta = zeros((qfM, MM * N));
θ2D_vartheta = zeros((qfM, MM * N));

for i in 0:MM * N-1
    #v odd that this is required. but otherwise the mod function removes some of the input???
    arr = 0:qfM
    idx = @. mod(arr - i * p, qfM) .+ 1

    idx = idx[1:end-1]
    #println(idx)

    #this will not work.
    r2D_vartheta[:, i+1] = r2D_alpha[idx, i+1]
    θ2D_vartheta[:, i+1] = θ2D_alpha[idx, i+1]
end
#%%

#just fourier transform of the above, doesn't add very much tbh.
#not super clear why we do this, maybe interpolating the fourier coeffs is easier?
scos_surf, ssin_surf = MID.QFM.rfft2D(r2D_vartheta, M, N);
tcos_surf, tsin_surf = MID.QFM.rfft2D(θ2D_vartheta, M, N);
#%%
contourf(tsin_surf)

#%%
contourf(r2D_alpha)
contourf(r2D_vartheta)
contourf(θ2D_alpha')
contourf(θ2D_vartheta')
#%%

#so it looks like the field lines are basically just extended to the new domain. 
#i.e. we are going from 2π to 2qπ. Unsure why we are doing this tbh.
display(size(r2D_alpha))
contourf(r2D_alpha)
contourf(r)
contourf(θ2D_alpha)
contourf(θ)

#%%
display(size(r))
p = plot()
for i in 1:32
    #plot!(r[i, :])
    plot!(θ[i, :])
end
display(p)
#%%
plot(θ)
contourf(θ')
plot(rcosarr[2, :])
plot(rsinarr)
plot(θcosarr)
plot(θsinarr)

#%%
ρ_surf = scos_surf[1, 1] #surface label.
surfaces = MID.QFM.QFMSurfaceT[]
push!(surfaces, MID.QFM.QFMSurfaceT(ρ_surf, scos_surf, tsin_surf, ssin_surf, tcos_surf))
plot_surfs(surfaces)

