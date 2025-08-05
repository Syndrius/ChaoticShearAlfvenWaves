
#actual integration with MID
using MID
using Plots
#%%

flr = MID.Structures.FLRT(δ=1e-8)
geo = init_geo(R0=1.0)
prob = MID.Helmholtz.TestProblemT(flr=flr, geo=geo)

Nx = 5
Ny = 5
Nz = 5
#this could cause problemos.
xgrid = init_grid(type=:rf, N = Nx, gp=3)
#so with pf = 3, we are very accuratly able to get the m=3 eigenvalues
#as we would hope,
#ideally, we can actually study this behaviour
#in order to make sure it is working as intended.
#but this looks pretty good.
#similar is noted for ζ.
#think we need to pick a relatively random assortment of evals 
#and see if we can hone in on them 
#perhaps consider the first m and second n or something?
#for each eval we should be able to pretty easily determine the zero by just counting the number of times the functions flips sign.
#guess we actually need to consider how important this all is?
#maybe better of just picking a few evals, 
#then looking at the convergence as res increases in parallel?
#This should show us that the grid is split properly
#but also can check pf by just doing this multiple times and hopefully showing that the convergence is faster for the chosen pf.
#last thing we really ought to check is different grid widths, eg run with quadratic etc.
ygrid = init_grid(type=:af, N = Ny, gp=3)#, pf=3)
zgrid = init_grid(type=:af, N = Nz, gp=3)#, pf=2)
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
evals, ϕ, ϕft = MID.Solve.compute_spectrum(prob=prob, grids=grids, solver=solver, deriv=true);
evals.ω[35:50] .^ 2
#%%
nbcs = length(unique(MID.Indexing.compute_boundary_inds(grids)))
ems = sort(real.(evals.ω .^2))
sort(real.(evals.ω))[nbcs+1:nbcs+20] .^2
nbcs = length(unique(MID.Indexing.compute_boundary_inds(grids)))
ems[nbcs+1:end]
#%%
xgp = MID.inst_grid(xgrid);
ygp = MID.inst_grid(ygrid);
zgp = MID.inst_grid(zgrid);
ygp = LinRange(0, 2π, MID.ifft_size(ygrid)+1)[1:end-1]
zgp = LinRange(0, 2π, MID.ifft_size(zgrid)+1)[1:end-1]

#for zero crossings we may want to focus on one of the fss cases!
ev = nbcs+37
display(ems[ev])
plot(xgp, real.(ϕft[ev, :, 1, 3, 1]))
plot!(xgp, real.(ϕft[ev, :, 1, 3, 5]))
#think to compare with analytical solution we might need to automatically find the phase.
contourf(ygp, xgp, real.(ϕ[ev, :, :, 2, 1]))
contourf(zgp, xgp, real.(ϕ[ev, :, 2, :, 1]))
contourf(zgp, ygp, real.(ϕ[ev, 2, :, :, 1]))

#%%
#finding number of zeros!
#could work, sign of m, n will be annoying
#and we might still have problemos with larger n.
nflip = 0
for i in 2:Nx-2 #ignore start and finish as we know they are zeros
    if sign(real(ϕft[ev, i, 1, 3, 1])) ≠ sign(real(ϕft[ev, i+1, 1, 3, 1])) #obvs will need the m, n combo we are using! -> annoying for sign of m, n
        nflip += 1
    end
end
display(nflip)
#%%
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
ktot, code = anal_evals(1:500, -100:100, -100:100);
(sort(real.(evals.ω))[nbcs+1:nbcs+20]) .^2
sort(real.(evals.ω)) .^2
sort(ktot)[1:20]
sort(real.(evals.ω))[1:20] .^2

length(ktot[ ktot .< 100])


#%%
#I think we just need to pick our ~10 favourite modes.
#this will probably require ~10 in each dim. so parallel.
#ideally we would be able to just find that these are the 10th, 15th etc eigenmode
#but I don't think that will be true,
#at least for low res.
#we may have to hunt through the solutions.
#specifically the m=0,n=0 case seems a bit odd, but perhaps not now we aren't doing the extra deriv shite.
#think the early ones will be indiciable, but no the larger ones.
#guess we can just go from the analytical ones and above
#to find the efunc that has the same m, n and number of zeros.
mlist = [1, 3, 5]
nlist = [1, 3, 5]
zl = [1, 3, 5]
mlist = [1, 3]
nlist = [1, 3]
zl = [1, 3]
#these should all be resolvable even in the small case.
mlist = [1, 2]
nlist = [1, 2]
zl = [1, 2]

#actually not baahd, just obvs need to factor in which zero is the solution!
#can either actually check the efunc, or maybe just group the evals somehow? Think that will be numerically difficult.
for m in mlist, n in nlist, zr in zl
    display("True sol")
    display(anal_evals(zr:zr, m:m, n:n))
    mat = []
    for i in nbcs+1:nbcs+200
        if (m, n) == (abs(evals.modelabs[i][1]), abs(evals.modelabs[i][2]))
            push!(mat, i)
            display(evals.ω[i] ^2)
        end
    end
    zrs = []
    zr_count = 1
    push!(zrs, mat[1])
    for i in mat
        #because the n label is not reliable, this does not really work!
        #might actually have to check how many times the sign changes.
        #this will naturally become less of a problemo with higher res
        #or with fss.
        if real(evals.ω[i] - evals.ω[zrs[zr_count]]) > 0.5 #defs a new zero
            push!(zrs, i)
            zr_count += 1
        end
    end
    display("Found")
    display(evals.ω[zrs[zr]]^2)

end

#this seems to actually be making some sense now!
real.(evals.ω)[nbcs+30:nbcs+50] .^ 2

evals.modelabs[nbcs+1:nbcs+20]




