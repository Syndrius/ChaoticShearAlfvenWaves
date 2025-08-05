
#first actual convergence test with the real case.
#don't think we can realisticly do this, as higher grids results in more
#k vals, so I think we are better off just looking at a tae frequency as the grid resolution goes up.
using MID
using Plots; plotlyjs()
using MIDViz
using FFTW
include("Solution.jl")
#%%

function fss_convergence(Nfs, Nss, mtarg, ntarg, zrtarg)

    #assume Nfs and Nsf are hte same length, just the size of each grid.

    flr = MID.Structures.FLRT(δ=1e-8)
    geo = init_geo(R0=1.0)
    prob = MID.Helmholtz.TestProblemT(flr=flr, geo=geo)
    an_ev, code = anal_evals(mtarg, ntarg, zrtarg)
    solver = init_solver(nev=20, targets=sqrt.(an_ev), prob=prob)
    eval_diff = zeros(length(Nfs), length(an_ev))
    efunc_diff = zeros(length(Nfs), length(an_ev))

    Nint = 2 * Nfs[end]
    int_grid = LinRange(0, 1, Nint)

    an_sol = zeros(length(an_ev), Nint)
    intphi = zeros(ComplexF64, length(an_ev), Nint) #needs to be complex unfort
    
    #gets the radial part of each eigenvalue.
    #comparing the angular parts is more effort than worth.
    for i in 1:length(an_ev)
        an_sol[i, :] = rad_anal_sol(code[i][1], code[i][2], int_grid)
    end

    for i in 1:length(Nfs)
        xgrid = init_grid(type=:rf, N = Nfs[i])#, gp=3)
        ygrid = init_grid(type=:as, start=1, N=Nss[i])
        zgrid = init_grid(type=:as, start=1, N=Nss[i])
        
        grids = init_grids(xgrid, ygrid, zgrid)
        evals, ϕ, ϕft = MID.Solve.compute_spectrum(prob=prob, grids=grids, solver=solver, deriv=true);
        #don't need this with slice solver
        #nbcs = length(unique(MID.Indexing.compute_boundary_inds(grids)))

        eval_true = []
        mn_true = []
        for zr in zrtarg, m in mtarg, n in ntarg
            mat = []
            mn = []
            for j in 1:length(evals.ω)
                if (m, n) == (abs(evals.modelabs[j][1]), abs(evals.modelabs[j][2]))
                    push!(mat, j)
                    push!(mn, evals.modelabs[j])
                end
            end
            for (j, ind) in enumerate(mat)
                #only works with m, n starting from 1.
                if zero_crossings(ϕft[ind, :, m, n, 1], Nfs[i]) == zr

                    push!(eval_true, ind)
                    push!(mn_true, mn[j])
                    break
                end
            end
        end
        eval_diff[i, :] = abs.(sort(real.(evals.ω[eval_true])) .^2 .- sort(an_ev))
        #now lets see if we can compare the efuncs as well.
        #probably easiest to just consider the radial part.
        #for this we will interpolate the solution, as the non-interpolated parts should match exactly.

        xgp = MID.inst_grid(xgrid)
        for k in 1:length(eval_true), j in 1:Nint
            #will only work if m,n are positive incrementing from 1.
            intphi[k, j] = MID.Mapping.hermite_interpolation(int_grid[j], ϕft[eval_true[k], :, mn_true[k][1], mn_true[k][2], :], xgp)
        end

        for k in 1:length(eval_true)
            am1 = argmax(abs.(an_sol[k, :]))
            am2 = argmax(abs.(real.(intphi[k, :])))
            sf = an_sol[k, am1] / real(intphi[k, am2])
            efunc_diff[i, k] = sum(abs.(an_sol[k, :] .- sf .* real.(intphi[k, :])))
        end

    end
    return eval_diff, efunc_diff
end

function ffs_convergence(Nfs, Nss, mtarg, ntarg, zrtarg)

    #assume Nfs and Nsf are hte same length, just the size of each grid.

    flr = MID.Structures.FLRT(δ=1e-8)
    geo = init_geo(R0=1.0)
    prob = MID.Helmholtz.TestProblemT(flr=flr, geo=geo)
    an_ev, _ = anal_evals(mtarg, ntarg, zrtarg)
    solver = init_solver(nev=30, targets=sqrt.(an_ev), prob=prob)
    eval_diff = zeros(length(Nfs), length(an_ev))
    efunc_diff = zeros(length(Nfs), length(an_ev))

    Nint = 2 * Nfs[end]
    int_grid = LinRange(0, 1, Nint)

    an_sol = zeros(length(an_ev), Nint)
    intphi = zeros(ComplexF64, length(an_ev), Nint, Nint) #needs to be complex unfort
    
    #gets the radial part of each eigenvalue.
    #comparing the angular parts is more effort than worth.
    for i in 1:length(an_ev)
        an_sol[i, :] = rad_anal_sol(code[i][1], code[i][2], int_grid)
    end

    for i in 1:length(Nfs)
        xgrid = init_grid(type=:rf, N = Nfs[i])#, gp=3)
        ygrid = init_grid(type=:af, N=Nfs[i])
        zgrid = init_grid(type=:as, start=1, N=Nss[i])
        
        grids = init_grids(xgrid, ygrid, zgrid)
        evals, ϕ, ϕft = MID.Solve.compute_spectrum(prob=prob, grids=grids, solver=solver, deriv=true);
        #nbcs = length(unique(MID.Indexing.compute_boundary_inds(grids)))

        eval_true = []
        mn_true = []
        for zr in zrtarg, m in mtarg, n in ntarg
            mat = []
            mn = []
            for j in 1:length(evals.ω)
                if abs((abs(evals.ω[j]) - 1)) < 0.2 #bc solution, so we ignore.
                    continue
                end
                if (m, n) == (abs(evals.modelabs[j][1]), abs(evals.modelabs[j][2]))
                    push!(mat, j)
                    push!(mn, evals.modelabs[j])
                end
            end
            for (j, ind) in enumerate(mat)
                #only works with m, n starting from 1.
                #this will get more convaluted with more fem
                #this needs to be done better, i.e. need to pick m, n better.
                #m starts from zero
                if zero_crossings(ϕft[ind, :, m+1, n, 1], Nfs[i]) == zr

                    push!(eval_true, ind)
                    push!(mn_true, mn[j])
                    break
                end
            end
        end
        eval_diff[i, :] = abs.(sort(real.(evals.ω[eval_true])) .^2 .- sort(an_ev))
        xgp = MID.inst_grid(xgrid)
        ygp = MID.inst_grid(ygrid)
        for k in 1:length(eval_true), j in 1:Nint, l in 1:Nint
            #will only work if m,n are positive incrementing from 1.
            #now m, n will be extra cooked.
            #ah yes, we have -1, need some kind of ind_to_mode map.
            #no longer actually needed as we need to interpolate in 2d.
            #mtrue = mode_to_ind(mn_true[k][1], ygrid)
            #mn_true for the last bit doesn't make any sense!
            intphi[k, j, l] = MID.Mapping.hermite_interpolation(int_grid[j], int_grid[l], ϕ[eval_true[k], :, :, mn_true[k][2], :], xgp, ygp)
        end

        intphi = fft(intphi, [3])

        for k in 1:length(eval_true)
            #ntrue = mn_true[k][2] #don't think this is needed, just because we are essentially picking a random part.
            am1 = argmax(abs.(an_sol[k, :]))
            am2 = argmax(abs.(real.(intphi[k, :, 1])))
            sf = an_sol[k, am1] / real(intphi[k, am2, 1])
            efunc_diff[i, k] = sum(abs.(an_sol[k, :] .- sf .* real.(intphi[k, :, 1])))
        end

    end
    return eval_diff, efunc_diff
end

#%%
function fff_convergence(Nfs, mtarg, ntarg, zrtarg)

    #assume Nfs and Nsf are hte same length, just the size of each grid.

    flr = MID.Structures.FLRT(δ=1e-8)
    geo = init_geo(R0=1.0)
    prob = MID.Helmholtz.TestProblemT(flr=flr, geo=geo)
    an_ev, _ = anal_evals(mtarg, ntarg, zrtarg)
    solver = init_solver(nev=30, targets=sqrt.(an_ev), prob=prob)
    eval_diff = zeros(length(Nfs), length(an_ev))
    efunc_diff = zeros(length(Nfs), length(an_ev))

    Nint = 2 * Nfs[end]
    int_grid = LinRange(0, 1, Nint)

    an_sol = zeros(length(an_ev), Nint)
    intphi = zeros(ComplexF64, length(an_ev), Nint, Nint, Nint) #needs to be complex unfort
    
    #gets the radial part of each eigenvalue.
    #comparing the angular parts is more effort than worth.
    for i in 1:length(an_ev)
        an_sol[i, :] = rad_anal_sol(code[i][1], code[i][2], int_grid)
    end

    for i in 1:length(Nfs)
        xgrid = init_grid(type=:rf, N = Nfs[i])#, gp=3)
        ygrid = init_grid(type=:af, N=Nfs[i])
        zgrid = init_grid(type=:af, N=Nfs[i])
        
        grids = init_grids(xgrid, ygrid, zgrid)
        evals, ϕ, ϕft = MID.Solve.compute_spectrum(prob=prob, grids=grids, solver=solver, deriv=true);
        #nbcs = length(unique(MID.Indexing.compute_boundary_inds(grids)))

        eval_true = []
        mn_true = []
        for zr in zrtarg, m in mtarg, n in ntarg
            display((zr, m, n))
            mat = []
            mn = []
            for j in 1:length(evals.ω)
                if abs((abs(evals.ω[j]) - 1)) < 0.2 #bc solution, so we ignore.
                    continue
                end
                if (m, n) == (abs(evals.modelabs[j][1]), abs(evals.modelabs[j][2]))
                    push!(mat, j)
                    push!(mn, evals.modelabs[j])
                end
            end
            for (j, ind) in enumerate(mat)
                #only works with m, n starting from 1.
                #this will get more convaluted with more fem
                #this needs to be done better, i.e. need to pick m, n better.
                #m starts from zero

                #mind = mode_to_ind(mn[j][1], ygrid)
                #nind = mode_to_ind(mn[j][2], zgrid)
                #display("here")
                #think the m+1 and n+1 is only valid if our analytical cases start from 1.
                #display(size(ϕft))
                #display((m+1, n+1))
                if zero_crossings(ϕft[ind, :, m+1, n+1, 1], Nfs[i]) == zr

                    push!(eval_true, ind)
                    push!(mn_true, mn[j])
                    #display("or here perhaps?")
                    break
                end
            end
        end
        #display("or here")
        eval_diff[i, :] = abs.(sort(real.(evals.ω[eval_true])) .^2 .- sort(an_ev))

        xgp = MID.inst_grid(xgrid)
        ygp = MID.inst_grid(ygrid)
        zgp = MID.inst_grid(zgrid)
        for k in 1:length(eval_true), j in 1:Nint, l in 1:Nint, p in 1:Nint
            #will only work if m,n are positive incrementing from 1.
            #now m, n will be extra cooked.
            #ah yes, we have -1, need some kind of ind_to_mode map.
            #no longer actually needed as we need to interpolate in 2d.
            #mtrue = mode_to_ind(mn_true[k][1], ygrid)
            #mn_true for the last bit doesn't make any sense!
            intphi[k, j, l, p] = MID.Mapping.hermite_interpolation(int_grid[j], int_grid[l], int_grid[p], ϕ[eval_true[k], :, :, :, :], xgp, ygp, zgp)
        end

        intphi = fft(intphi, [3, 4])

        for k in 1:length(eval_true)
            #ntrue = mn_true[k][2] #don't think this is needed, just because we are essentially picking a random part.
            mtrue, ntrue = mn_true[k]
            if mtrue < 0
                mind = Nint + mtrue + 1
            else
                mind = mtrue + 1
            end
            if ntrue < 0
                nind = Nint + ntrue + 1
            else
                nind = ntrue + 1
            end
            am1 = argmax(abs.(an_sol[k, :]))
            am2 = argmax(abs.(real.(intphi[k, :, mind, nind])))
            sf = an_sol[k, am1] / real(intphi[k, am2, mind, nind])
            efunc_diff[i, k] = sum(abs.(an_sol[k, :] .- sf .* real.(intphi[k, :, mind, nind])))
        end

    end
    return eval_diff, efunc_diff

end
#%%
#think this is working!
#this should defs be in MID.
function mode_to_ind(m::Int64, grid::MID.Structures.AngFEMGridDataT)
    #going to assume we don't have a pf for now!
    
    if m < 0
        ind = grid.N + m + 1
    else
        ind = m + 1 #account for zero?
    end
    return ind
end
tgrid = init_grid(type=:af, N=5)
mode_to_ind(2, tgrid)
#%%
#need to compare the actual eigen functions as well.
#and be able to do this with weird grids and the pf,
#although basic testing shows pf working, concern is m vs -m in our actual code
#which is ideantical in this case.

#looks to be working, under a few constraints.
#actually choosing the target evals will be v difficult, especially for fem.
#actually seems to be working with the efunc!
dval, dfunc = fss_convergence([10, 15], [2, 3], [1, 2], [1, 2], [1, 2])
#now working, looks good, ideally we can compare the error for different m's or something.
#efuncs look ok, bit skeptical that correct modes are actually being chosen, but it doesn't really matter.
dval, dfunc = ffs_convergence([10, 15], [2, 3], [1, 2], [1, 2], [1, 2])
#also looks to be working
#this is slow af, need to parallise for this to be practical!
#actually also looks good, clear problemos when a specific mode cannot be resolved.
dval, dfunc = fff_convergence([6, 8], [1, 2], [1, 2], [1, 2])
#%%
p = plot()
for i in 1:2
    scatter!(dval[i, :])
end
display(p)
#%%
p = plot()
for i in 1:2
    scatter!(dfunc[i, :])
end
display(p)
#%%
evals.ω .^ 2  
mlist = [1, 2]
nlist = [1, 2]
zl = [1, 2]

#actually not baahd, just obvs need to factor in which zero is the solution!
#can either actually check the efunc, or maybe just group the evals somehow? Think that will be numerically difficult.
eval_true = []
for m in mlist, n in nlist, zr in zl
    display("True sol")
    display(anal_evals(zr:zr, m:m, n:n))
    mat = []
    mn = []
    for i in nbcs+1:length(evals.ω)
        if (m, n) == (abs(evals.modelabs[i][1]), abs(evals.modelabs[i][2]))
            push!(mat, i)
            push!(mn, evals.modelabs[i])
        end
    end
    zrs = []
    zr_count = 1
    push!(zrs, mat[1])
    for i in mat
        #only works with m, n starting from 1.
        if zero_crossings(ϕft[i, :, m, n, 1], Nx) == zr

            push!(eval_true, i)
            break
        end
    end

        #because the n label is not reliable, this does not really work!
        #might actually have to check how many times the sign changes.
        #this will naturally become less of a problemo with higher res
        #or with fss.
        #if real(evals.ω[i] - evals.ω[zrs[zr_count]]) > 0.5 #defs a new zero
        #    push!(zrs, i)
        #    zr_count += 1
        #end
    display("Found")
    display(evals.ω[eval_true[end]]^2)

end
#%%
#coolio, we are seeing convergence already!
#naturally we need to run this with larger choices for m, n and zr.
anev, _ = anal_evals(mlist, nlist, zl)

sort(abs.(evals.ω[eval_true] .^2)) .- sort(anev)
#%%
Nfs = [10]
Nss = [2]
mtarg = [1, 2]
ntarg = [1, 2]
zrtarg = [1, 2]
flr = MID.Structures.FLRT(δ=1e-8)
geo = init_geo(R0=1.0)
prob = MID.Helmholtz.TestProblemT(flr=flr, geo=geo)
an_ev, code = anal_evals(mtarg, ntarg, zrtarg)
solver = init_solver(nev=20, targets=sqrt.(an_ev), prob=prob)
eval_diff = zeros(length(Nfs), length(an_ev))

xgrid = init_grid(type=:rf, N = 10)
ygrid = init_grid(type=:as, start=1, N=2)
ygrid = init_grid(type=:af, N = 10)
zgrid = init_grid(type=:as, start=1, N=2)

grids = init_grids(xgrid, ygrid, zgrid)
evals, ϕ, ϕft = MID.Solve.compute_spectrum(prob=prob, grids=grids, solver=solver, deriv=true);
#nbcs = length(unique(MID.Indexing.compute_boundary_inds(grids)))
Nint = 2 * Nfs[end]
int_grid = LinRange(0, 1, Nint)

an_sol = zeros(length(an_ev), Nint)

perm = sortperm(an_ev)
san_ev = an_ev[perm]
scode = code[perm]
#gets the radial part of each eigenvalue.
#comparing the angular parts is more effort than worth.
for i in 1:length(an_ev)
    #should be in sorted order!
    #an_sol[i, :] = rad_anal_sol(scode[i][1], scode[i][2], int_grid)
    an_sol[i, :] = rad_anal_sol(code[i][1], code[i][2], int_grid)
end

eval_true = []
#don't need to worry about sorting if we just loop thorugh this in the same order as the analytical ones are found.
for zr in zrtarg, m in mtarg, n in ntarg
    mat = []
    mn = []
    for j in 1:length(evals.ω)
        if abs((abs(evals.ω[j]) - 1)) < 0.2 #bc solution
            continue
        end
        if (m, n) == (abs(evals.modelabs[j][1]), abs(evals.modelabs[j][2]))
            push!(mat, j)
            push!(mn, evals.modelabs[j])
        end
    end
    for j in mat
        #if zero_crossings(ϕft[j, :, m, n, 1], 10) == zr
        if zero_crossings(ϕft[j, :, m+1, n, 1], 10) == zr
            display((m, m, zr))
            push!(eval_true, j)
            break
        end
    end
end
display(eval_true)
abs.(sort(real.(evals.ω[eval_true])) .^2 .- sort(an_ev))
#%%
am1 = argmax(abs.(an_sol[7, :]))
am2 = argmax(abs.(real.(intphi)))
sf = an_sol[7, am1] / real(intphi[am2])
#potential_plot(ϕft, grids, eval_true[7], label_max=0.5)
plot(xgp, real.(ϕft[eval_true[7], :, 3, 1, 1]))
plot!(int_grid, an_sol[7, :])
plot!(int_grid, sf .* real.(intphi[:, 1]))# sort of works!
#%%
#so 1d interpolation doesn't really work unfort.
#ffs case is actually awkward af, we probably need the original ϕ to make this work, i.e. ffs shape.
#guess we can just assume a sin and hope for the best, i.e. that we havent picked zero lol
intphi = zeros(ComplexF64, Nint, Nint)
ygp = MID.inst_grid(ygrid)
for i in 1:Nint, j in 1:Nint
    intphi[i, j] = MID.Mapping.hermite_interpolation(int_grid[i], int_grid[j], ϕ[eval_true[7], :, :, 1, :], xgp, ygp)
end

intphi = fft(intphi, [2])
#%%
#so 1d interpolation doesn't really work unfort.
intphi = zeros(ComplexF64, Nint)

for i in 1:Nint
    intphi[i] = MID.Mapping.hermite_interpolation(int_grid[i], ϕft[eval_true[7], :, 3, 1, :], xgp)
end
#%%
xgp = MID.inst_grid(xgrid)
plot(xgp, real.(ϕft[eval_true[2], :, 2, 1, 1])) 

eval_true
#direct example of modelabs being wrong!
evals.modelabs[2]
sort(real.(evals.ω .^2))
