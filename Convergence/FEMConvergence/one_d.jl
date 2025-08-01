
using MID
using Statistics
using Plots
#%%

N = 10
xgrid = LinRange(0, π, N)
#xgrid = LinRange(0, 1, N)
gp = 3
llam, lf = MID.Helmholtz.linear(N, gp)

clam, cf = MID.Helmholtz.cubic(N, gp)

clam
#%%
xgrid = LinRange(0, 1, N)
plot(xgrid, real.(cf[1:2:end, 3]))
plot!(xgrid, real.(cf[2:2:end, 3]))

#%%

kex = 2

#first two are cooked beyond anything, but the others are good for linear case
as = sin.(xgrid .* kex)
sf = as[3] / lf[3, kex+2]
plot(xgrid, sf .*lf[:, kex+2])
plot!(xgrid, sin.(xgrid .* kex))#*π))
#%%
#looks like amp is just not normalised, I guess any amplitude would satisfy the equaiton, hence they have
#to be normalised.
#unsure how the fek to do that...
#I guess we can do that with the analytical solution?
#I suppose we would need better boundary conditions for this to actually have a unique solution.
as = sin.(xgrid .* kex)
sf = as[3] / cf[5, kex+2]

plot(xgrid, sf .*cf[1:2:end, kex+2])
plot!(xgrid, sin.(xgrid .* kex))

#%%
#ok so proper convergence test based on analytical solution.

#guess we will just look at the first few k values?

Nlist = [10, 16, 32, 64, 128]

ktest = [3, 4, 5, 6, 7, 8]
ktrue = (ktest .- 2) .^ 2

lkerr = zeros(length(Nlist), length(ktest));
ckerr = zeros(length(Nlist), length(ktest));
lferr = zeros(length(Nlist), length(ktest));
cferr = zeros(length(Nlist), length(ktest));
gp = 4
cfs = []
for (i, N) in enumerate(Nlist)
    xgrid = LinRange(0, π, N)
    llam, lf = MID.Helmholtz.linear(N, gp)
    clam, cf = MID.Helmholtz.cubic(N, gp)
    
    lkerr[i, :] = abs.(llam[ktest] .- ktrue)
    ckerr[i, :] = abs.(clam[ktest] .- ktrue)

    for (j, k) in enumerate(ktest)
        ftrue = sin.(k .* xgrid)
        lsf = ftrue[3] / lf[3, 2+k]
        csf = ftrue[3] / cf[5, 2+k]
        lferr[i, j] = mean(abs.(ftrue .- lsf .* lf[:, 2+k]))
        cferr[i, j] = mean(abs.(ftrue .- csf .* cf[1:2:end, 2+k]))
        push!(cfs, csf .* cf[1:2:end, 2+k])
    end

end
#%%

p = plot()
for i in 1:length(Nlist)

    xgrid = LinRange(0, π, Nlist[i])

    ftrue = sin.(3*xgrid)
    plot!(xgrid, cfs[i*6])#[1:2:end, 5])
end
display(p)

#scatter(lkerr[:, 4])
#scatter!(ckerr[:, 4])
#%%
#we probably have to interpolate between the values for this to mean anything!
#but I guess this is saying that at the exact points the solution matches to 10^-13.
#partly because we have scaled it that way.
#I guess for 3d we can probably check the interpolated case.
p = plot()#ylimits=(0, 1))
for i in 1:length(ktest)
    #scatter!(log.(lkerr[:, i]))
    #scatter!(log.(ckerr[:, i]))
    #scatter!((lferr[2:end, i]))
    scatter!((cferr[1:end, i]))
end

display(p)


#%%
#integration is wrong yip fkn yee
#lets build up the single element K matrix.
using FastGaussQuadrature
ξ, wg = gausslegendre(gp)

S = MID.Helmholtz.linear_basis(ξ)

#this should just have a dx factor.
K = zeros(2, 2)
I = zeros(2, 2)

for i in 1:gp
    K[1, 1] += S.dH[1, i] * S.dH[1, i] * wg[i]
    K[1, 2] += S.dH[1, i] * S.dH[2, i] * wg[i]
    K[2, 1] += S.dH[2, i] * S.dH[1, i] * wg[i]
    K[2, 2] += S.dH[2, i] * S.dH[2, i] * wg[i]

    I[1, 1] += S.H[1, i] * S.H[1, i] * wg[i]
    I[1, 2] += S.H[1, i] * S.H[2, i] * wg[i]
    I[2, 1] += S.H[2, i] * S.H[1, i] * wg[i]
    I[2, 2] += S.H[2, i] * S.H[2, i] * wg[i]

end

#so integration is ok, we are somhow just doubling answer in second matrix?
display(K) #missing factor of 2, comes from jac / jac^2.
display(I) #missing factor of half, comes from jac.


