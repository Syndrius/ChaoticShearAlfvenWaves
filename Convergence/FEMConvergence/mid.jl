
#actual integration with MID
using MID
#%%

flr = MID.Structures.FLRT(δ=1e-8)
geo = init_geo(R0=1.0)
prob = MID.Helmholtz.TestProblemT(flr=flr, geo=geo)

Nx = 9
Ny = 5
Nz = 5
#this could cause problemos.
xgrid = init_grid(type=:rf, N = Nx, gp=3)
ygrid = init_grid(type=:af, N = Ny, gp=3)#, pf=2)
zgrid = init_grid(type=:af, N = Nz, gp=3)
#this seems to be find crossing zero though!
#ygrid = init_grid(type=:as, start=-1, N=3)
#looks like this doesn't work when n crosses zero.
#zgrid = init_grid(type=:as, start=-2, N=3)
#I think, because we have fourth order ζ term, n cannot be zero, as this cooks the equation.
#and yest, we seem to be able to find the n=0 eigenvalue sometimes.
#I guess specifically the matrices will be indeterminant for this case I think.
#zgrid = init_grid(type=:as, start=0, N=1)

grids = init_grids(xgrid, ygrid, zgrid)

solver = init_solver(target=sqrt(26), nev=100, prob=prob)
solver = init_solver(full_spectrum=true, prob=prob)

#%%

#this actually looks like it is working now...
#pretty hekin wild.
#the eval pattern will require a proper undertsanding of the solution.
#then we need to look at the ocnvergence as the grid size increases
#compare fss to fff -> may be a bit tricky as the modes are decoupled.
#make sure this works in v large scales in Petsc.
evals, ϕ, ϕft = MID.Solve.compute_spectrum(prob=prob, grids=grids, solver=solver);
#%%
nbcs = length(unique(MID.Indexing.compute_boundary_inds(grids)))
ems = sort(real.(evals.ω .^2))
sort(real.(evals.ω))[nbcs+1:nbcs+20] .^2
nbcs = length(unique(MID.Indexing.compute_boundary_inds(grids)))
ems[nbcs+1:end]
#%%
xgp = MID.inst_grid(xgrid);
ygp = MID.inst_grid(ygrid);
ygp = LinRange(0, 2π, MID.ifft_size(ygrid)+1)[1:end-1]
zgp = LinRange(0, 2π, MID.ifft_size(zgrid)+1)[1:end-1]

ev = nbcs+3
display(ems[ev])
plot(xgp, real.(ϕft[ev, :, 2, 1]))
#think to compare with analytical solution we might need to automatically find the phase.
contourf(ygp, xgp, real.(ϕ[ev, :, :, 2]))
contourf(zgp, xgp, real.(ϕ[ev, :, 2, :]))

#so ϕ goes 
#[efunc, rad, bessel function, n^2]
#eval is zero of bessel function + n^2
#n is unrealted to order of bessel function

#looks like mo=0 case is ok, Bessel function just gets distored to 0 at origina due to boundary condition.
#I guess that solution is probably some weird combo of J and Y. -> this will be impossible to compare for convergence
#so we won't worry about it!
ktrue = anal_evals(1:5, -2:2, 1:2) #.^ 2
display(ems[9:20])

ftrue = anal_sol(1, 1, 1, xgp, ygp, zgp);

#solution is clearly correct
#going to be difficult to compare though!
#we may need to compare the fft result
#i.e. just make sure we have the same harmonics
#and then compare the radial part?
contourf(zgp, xgp, ftrue[:, 2, :])
contourf(ygp, xgp, ftrue[:, :, 2])

plot(xgp, ftrue[:, 2, 2])
plot!(xgp, real.(ϕ[ev, :, 2, 2]))
#%%
#evals.modelabs
diff = zeros(length(evals.ω) - nbcs)
for i in nbcs+1:length(evals.ω)
    #so problems arise when the eigenvalue is correct
    #but corresponds to an eigenfunction with m > Ny/2 or whatever.
    #unsure how this should be dealt with.
    #don't really understand how it is able to find this eval tbh.
    #i guess cos(6x) looks identical to cos(3x) when the grid size is only 3. etc.
    m, n = evals.modelabs[i]
    #display((m, n))
    #not really sure how many evals there actually are for each case.
    #perhaps for a given (m, n) we just need to find the zero such that the eval is closet.
    #zr = Int64(ceil(count(j->(abs(j[1])==abs(m) && abs(j[2]==abs(n))), evals.modelabs[nbcs+1:i]) / 2))
    ktest = anal_evals(1:100, m:m, n:n)

    ktrue = ktest[argmin(abs.(ktest .- real(evals.ω[i]^2)))]
    #display(ktrue)
    diff[i-nbcs] = real(evals.ω[i]^2) - ktrue
    #display(evals.ω[i]^2)
    #display(ktrue)
end

scatter(diff[1:100])
#%%
argmax(diff)

evals.ω[nbcs+1]^2

evals.modelabs[nbcs+13]


contourf(ygp, xgp, real.(ϕ[nbcs+13, :, :, 2]))
contourf(zgp, xgp, real.(ϕ[nbcs+13, :, 2, :]))

plot(xgp, real.(ϕ[nbcs+13, :, 2, 2]))


#probably easiest to just compare th 100 smallest evals, excluding zero ofc.
#the ol -ve n and -ve m seem to be causing problemos...
#unsure how to reliably get the eigenvalues
#perhaps this won't be a problem once the grid is larger.
#this is perhaps highlighting a real problemo.
#think as the grid gets larger this is not a problem.
ktot = anal_evals(1:200, -50:50, -50:50);
(sort(real.(evals.ω))[nbcs+1:nbcs+20]) .^2
sort(real.(evals.ω)) .^2
sort(ktot)[1:20]
sort(real.(evals.ω))[1:20] .^2
