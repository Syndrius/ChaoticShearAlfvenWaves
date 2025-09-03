
#fkn stupid as name
#seeing if we can judge our crisp a mode is based on 
#the relative amplitude of the second largest harmonic
#the second derivative at the peak.
using MID
using MIDViz
using Plots#; plotlyjs()

#%%

#think we need to label the modes etc for each of these cases.

#terrible name but no idea what to call it.
struct HarmonicsT
    fwhm :: Array{Float64} #full width half maximum of dominant harmonic
    dd :: Array{Float64} # second derivative at the peak location of the dominant harmonic
    harm2 :: Array{Float64} # relative peak amplitude of the second largest harmonic
    modelabs :: Array{Tuple{Int64, Int64}} # same as evalsT, so we can plot nicely
    x1 :: Array{Float64} # same as evalsT, the radial location of the dominant harmonic for plotting.
    function HarmonicsT(N::Int64)
        new(zeros(N), zeros(N), zeros(N), Array{Tuple{Int64, Int64}}(undef, N), zeros(N))
    end
end

function crossings(f, N)

    am = argmax(abs.(real.(f)))
    nf = real.(f) / real(f[am]) .- 0.5

    nflip = 0
    locs = []
    for i in 1:N-1 #ignore start and finish as they are zero. -> now kept for real edge case where spike is at the edge or start
        if sign(real(nf[i])) ≠ sign(real(nf[i+1]))
            nflip += 1
            #returns the index of the point to the left of the crossing
            push!(locs, i)
        end
    end
    return nflip, locs 
end
#%%
#finds some information about the harmonic structure of the computed spectrum.
#this will be annoying af to do in general
#will be a bit like mapping, going to have to re-process everything...
function harmonic_info(ϕft, grids, evals)

    x1m = zeros(Int64, grids.x2.N, grids.x3.N)
    nevals = length(evals.ω)
    ht = HarmonicsT(nevals)

    x1grid = MID.inst_grid(grids.x1)

    intsize = 10
    intg1 = zeros(intsize)
    intg2 = zeros(intsize)
    ϕ1 = zeros(ComplexF64, intsize)
    ϕ2 = zeros(ComplexF64, intsize)

    ϕm = zeros(grids.x2.N, grids.x3.N)
    nvals = size(ϕft)[1]
    ht.modelabs .= evals.modelabs

    for i in 1:nvals

        for j in 1:grids.x2.N, k in 1:grids.x3.N
            #note we do not need the derivs for this case~
            x1m[j, k] = argmax(abs.(real.(ϕft[i, :, j, k, 1])))
            ϕm[j, k] = abs.(real.(ϕft[i, x1m[j, k], j, k, 1]))
        end
        max_mode = argmax(ϕm)
        max_val = ϕm[max_mode]
        ϕm[max_mode] = 0.0 #so we can find the second
        max2 = argmax(ϕm)

        ht.harm2[i] = ϕm[max2] / max_val
        ht.x1[i] = x1grid[x1m[max_mode]]

        #computes the second derivative at the peak.
        #think this will be completly useless tbh!
        ht.dd[i] = real(MID.Mapping.hermite_interpolation(ht.x1[i], ϕft[i, :, max_mode, :], x1grid, 2))

        nflip, locs = crossings(ϕft[i, :, max_mode, 1], grids.x1.N)
        if length(locs) > 2
            locs = [locs[1], locs[2]]
        end

        intg1 .= LinRange(x1grid[locs[1]], x1grid[locs[1]+1], intsize)
        intg2 .= LinRange(x1grid[locs[2]], x1grid[locs[2]+1], intsize)

        for i in 1:intsize
            ϕ1[i] = MID.Mapping.hermite_interpolation(intg1[i], ϕft[i, :, max_mode, :], x1grid)
            ϕ2[i] = MID.Mapping.hermite_interpolation(intg2[i], ϕft[i, :, max_mode, :], x1grid)
        end

        ind1 = find_ind(ϕ1, 0.5)
        ind2 = find_ind(ϕ2, 0.5)
        ht.fwhm[i] = intg2[ind2] - intg1[ind1]
    end

    return ht
end

ht = harmonic_info(ϕft, grids, evals)
scatter(ht.x1, ht.fwhm)
scatter(ht.x1, ht.harm2)
scatter(ht.x1, ht.dd)
#%%


#conside a toroidal case, as we should see the second harmonic grow in size as a function of r.

geo = init_geo(R0=4.0)
prob = init_problem(q=fu_dam_q, geo=geo); 
#%%
#then create the grids
Nr = 100;
rgrid = init_grid(type=:rf, N=Nr)
θgrid = init_grid(type=:as, N = 3, start = 1)
ζgrid = init_grid(type=:as, N = 2, start = -2)
grids = init_grids(rgrid, θgrid, ζgrid);

#%%
#then define the solver
solver = init_solver(prob=prob, targets=[0.2, 0.3, 0.4], nev=80)
#%%

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver, deriv=true);
#%%
continuum_plot(evals)

ind = find_ind(evals, 0.2884)
harmonic_plot(ϕft, grids, ind)
#%%
#finally we want to do some kind of fwhm approx
#need to find how many times the main harmonic crosses 0.5 (after normalising!)
#then find the left most crossing and right most
#get closest two points
#interpolated between them
#then find points closest to 0.5
#and take the difference.

function crossings(f, N)

    am = argmax(abs.(real.(f)))
    nf = real.(f) / real(f[am]) .- 0.5

    nflip = 0
    locs = []
    for i in 1:N-1 #ignore start and finish as they are zero. -> now kept for real edge case where spike is at the edge or start
        if sign(real(nf[i])) ≠ sign(real(nf[i+1]))
            nflip += 1
            #returns the index of the point to the left of the crossing
            push!(locs, i)
        end
    end
    return nflip, locs 
end

ev = 213
#first get the radial profile of the dominant harmonic

x1m = zeros(Int64, grids.x2.N, grids.x3.N)
ϕm = zeros(grids.x2.N, grids.x3.N)

for j in 1:grids.x2.N, k in 1:grids.x3.N

    x1m[j, k] = argmax(abs.(real.(ϕft[ev, :, j, k, 1])))
    ϕm[j, k] = abs.(real.(ϕft[ev, x1m[j, k], j, k, 1]))
end

max_mode = argmax(ϕm)

plot(xgrid, real.(ϕft[ev, :, max_mode, 1]))
plot(xgrid, real.(ϕft[ev, :, max_mode, 1]) / maximum(abs.(real.(ϕft[ev, :, max_mode, 1]))) .- 0.5)
harmonic_plot(ϕft, grids, ev)
nflip, locs = crossings(ϕft[ev, :, max_mode, 1], grids.x1.N)

display(xgrid[locs])
#%%
#now we interpolate between them
#this may not even be necesary really, with fine grid we are probably pretty close to 0.5 but whatever.
if length(locs) > 2
    #get the largest distance,
    #this will probably have flaws.
    locs = [locs[1], locs[end]]
end
intg1 = LinRange(xgrid[locs[1]], xgrid[locs[1]+1], 20)
intg2 = LinRange(xgrid[locs[2]], xgrid[locs[2]+1], 20)

ϕ1 = zeros(ComplexF64, 20)
ϕ2 = zeros(ComplexF64, 20)
for i in 1:20
    ϕ1[i] = MID.Mapping.hermite_interpolation(intg1[i], ϕft[ev, :, max_mode, :], xgrid)
    ϕ2[i] = MID.Mapping.hermite_interpolation(intg2[i], ϕft[ev, :, max_mode, :], xgrid)
end

ind1 = find_ind(ϕ1, 0.5)
ind2 = find_ind(ϕ2, 0.5)
fwhm = intg2[ind2] - intg1[ind1]
#%%
function compute_fwhm(ϕft, grids)

    nevals = size(ϕft)[1]
    fwhm = zeros(nevals)
    xgrid = MID.inst_grid(grids.x1)

    intsize = 10
    intg1 = zeros(intsize)
    intg2 = zeros(intsize)
    ϕ1 = zeros(ComplexF64, intsize)
    ϕ2 = zeros(ComplexF64, intsize)
    xloc = zeros(Int64, nevals)

    for i in 1:nevals
        max_mode = argmax(abs.(real.(ϕft[i, :, :, :, 1])))

        xloc[i] = max_mode[1]
        max_mode = CartesianIndex((max_mode[2], max_mode[3]))

        nflip, locs = crossings(ϕft[i, :, max_mode, 1], grids.x1.N)
        #display(i)
        #display(locs)

        #this probably hasn't been tested properly yet!
        if length(locs) > 2
            locs = [locs[1], locs[2]]
        end

        intg1 .= LinRange(xgrid[locs[1]], xgrid[locs[1]+1], intsize)
        intg2 .= LinRange(xgrid[locs[2]], xgrid[locs[2]+1], intsize)

        for i in 1:intsize
            ϕ1[i] = MID.Mapping.hermite_interpolation(intg1[i], ϕft[ev, :, max_mode, :], xgrid)
            ϕ2[i] = MID.Mapping.hermite_interpolation(intg2[i], ϕft[ev, :, max_mode, :], xgrid)
        end

        ind1 = find_ind(ϕ1, 0.5)
        ind2 = find_ind(ϕ2, 0.5)
        fwhm[i] = intg2[ind2] - intg1[ind1]
    end
    return fwhm, xloc
end


fwhm, xlocs = compute_fwhm(ϕft, grids)

scatter(xgrid[xlocs], fwhm)

#%%
#now consider the second derivative at the peak.
#this may not make heaps of sense in the fss case.

#first make a super fine grid near the peak.
#this will be tricky in 3d!

ev = 92

evals.x1[ev]

xgrid = MID.inst_grid(rgrid)

peak_ind = find_ind(xgrid, evals.x1[ev])

lp = xgrid[peak_ind - 1]
rp = xgrid[peak_ind - 1]

int_grid = LinRange(lp, rp, 10) #unsure how much we need here!
#do we even need this? can we not just interpolate the second deriv at the actual point?

yz = argmax(abs.(real.(ϕft[ev, peak_ind, :, :, 1])))

display(ϕft[ev, peak_ind, yz, 2])

display(ϕft)


#seems reasonable, might need to double check the old mapping function.
#but it does look like as a first pass we can do this
#does giv pretty absurd results for continuum modes that are very spiky at the peak
#maybe this will make more sense with finer grid.
MID.Mapping.hermite_interpolation(xgrid[peak_ind], ϕft[ev, :, yz, :], xgrid, 2)

ϕft[ev, peak_ind, yz, 2]


harmonic_plot(ϕft, grids, ev)


#%%
x1m = zeros(Int64, grids.x2.N, grids.x3.N)
ϕm = zeros(grids.x2.N, grids.x3.N)

for j in 1:grids.x2.N, k in 1:grids.x3.N

    x1m[j, k] = argmax(abs.(real.(ϕft[92, :, j, k, 1])))
    ϕm[j, k] = abs.(real.(ϕft[92, x1m[j, k], j, k, 1]))
end

max_mode = argmax(ϕm)

maxval = ϕm[max_mode]
ϕm[max_mode] = 0.0

max2 = argmax(ϕm)

ϕm[max2] / maxval

#%%
#gets the relative amplitude of the second largest harmonic, not particularly helpful for generic toroidal case
#but perhaps for qfm chaos.
function second_harmonic(ϕft, grids)
    x1m = zeros(Int64, grids.x2.N, grids.x3.N)
    ϕm = zeros(grids.x2.N, grids.x3.N)
    nvals = size(ϕft)[1]

    harm2 = zeros(nvals)
    x1loc = zeros(nvals) #location of the main harmonic

    x1grid = MID.inst_grid(grids.x1)

    for i in 1:nvals

        for j in 1:grids.x2.N, k in 1:grids.x3.N
            #note we do not need the derivs for this case~
            x1m[j, k] = argmax(abs.(real.(ϕft[i, :, j, k, 1])))
            ϕm[j, k] = abs.(real.(ϕft[i, x1m[j, k], j, k, 1]))
        end
        max_mode = argmax(ϕm)
        max_val = ϕm[max_mode]
        ϕm[max_mode] = 0.0 #so we can find the second
        max2 = argmax(ϕm)

        harm2[i] = ϕm[max2] / max_val
        if harm2[i] > 1.0
            display(i)
        end
        x1loc[i] = x1grid[x1m[max_mode]]
    end
    return harm2, x1loc
end

harm2, rloc = second_harmonic(ϕft, grids)


#%%
#think this sort of makes sense now, pretty fkn stupid though, as the largest values are simply the avoided crossings!
#maybe this will offer some insight for the qfm case. unlikely though!
scatter(rloc, harm2)


