#see if we can determine a good qfm configuration, using only a single toroidal mode number
using MID
using MIDViz
using Plots; plotlyjs()
using Plots; gr()
#%%

function q_prof(r::Float64)
    #we are going to stick with the q-prof
    #n=-3 harmonics for perturbations
    #if we need to look at tae, we will need to change the density profile
    
    a = 1.05
    b = 0.5
    c = 9.0
    d = 0.0
    return a + b*r^2 + c*r^6, 2 * b * r + 6*c*r^5
    #return a + b*r^2 + c*r^10 + d*r^4, 2*b*r + 10*c*r^9 + 4*d*r^3
    #return a + b*r^6, 6 * b * r^5
    #return a + b*r^4, 4 * b * r^3
end
function dens_prof(r::Float64)

    return 1.0 #- tanh((r-0.85)/0.15)
end
#%%

rgrid = init_grid(type=:rc, N = 250)
θgrid = init_grid(type=:as, start=0, N=10)
ζgrid = init_grid(type=:as, start=-4, N=4)

cont_grids = init_grids(rgrid, θgrid, ζgrid)

#%%

geo = init_geo(R0=4.0)
prob= init_problem(geo=geo, q=q_prof, dens=dens_prof)
#%%

evals_cont = compute_continuum(prob, cont_grids);

#%%

continuum_plot(evals_cont, cont_grids)

#%%
rgrid = init_grid(type=:rf, N=60)
θgrid = init_grid(type=:as, start=0, N=9)
ζgrid = init_grid(type=:as, start=-5, N=5)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%

solver = init_solver(nev = 200, targets=[0.25, 0.30, 0.35], prob=prob)

#%%
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);
#%%

continuum_plot(evals)
continuum_plot(evals, n=-2)

ind  = find_ind(evals, 0.2466243)

potential_plot(ϕft, grids, ind, label_max=0.5)

#%%
k = 0.0015
isl1 = init_island(m0=7, n0=-3, A=k/7)
isl2 = init_island(m0=8, n0=-3, A=k/8)
isls = [isl1, isl2]

isl_prob = init_problem(q=q_prof, dens=dens_prof, geo=geo, isls=isls)

#%%

Ntraj = 80
rlist = collect(LinRange(0.6, 0.85, Ntraj));
Nlaps = 500

poincare_plot(isl_prob, Nlaps, Ntraj, rlist);

#%%
#see if we can plot the edge of the gaps, this will help us make the tae gap 'straight'
function get_gaps(q, dens, n)

    better_q(r) = q(r)[1]

    sols = []
    freqs = []
    for m in 1:10
        om1(r) = (m/better_q(r) + n) / sqrt(dens(r))
        om2(r) = ((m+1)/better_q(r) + n) / sqrt(dens(r))

        #om1(r) = (m/better_q(r) + n) 
        #om2(r) = ((m+1)/better_q(r) + n)

        #display(om1(0.3))
        #display(om2(0.3))
        diff(r) = om1(r) + om2(r)
        sol = 0.0
        try
            sol = find_zero(diff, 0.5)
        catch

        end
        display(sol)
        push!(sols, sol)
        push!(freqs, om1(sol))
    end
    return sols, freqs
end 
#%%
function q_prof(r)
    a = 1.05
    b = 0.4
    c = 5.0
    d = 6
    return a + b*r^2 + c*r^d, 2*b*r + d*c*r^(d-1) 
    #return a + b*r^6, 6 * b * r^5
    #return a + b*r^4, 4 * b * r^3
    #return a + b*r^2, 2 * b * r
end
function dens_prof(r::Float64)
    a = 1.05
    b = 0.4
    c = 5.0
    d = 6

    if r < 0.6
        return 1.0
    end
    return 1/(-0.867647*r + 0.937941)^2
end

sols1, freqs1 = get_gaps(q_prof, dens_prof, -1)
sols2, freqs2 = get_gaps(q_prof, dens_prof, -2)
sols3, freqs3 = get_gaps(q_prof, dens_prof, -3)
sols4, freqs4 = get_gaps(q_prof, dens_prof, -4)
#%%

plot(sols1, abs.(freqs1), ylimits=(0.0, 1.5), xlimits=(0.0, 1.0))
plot!(sols2, abs.(freqs2))
#plot!(sols3, abs.(freqs3))
#plot!(sols4, abs.(freqs4))

#%%
rvals = LinRange(0, 1, 100)
omvals = [2 / q_prof(r)[1] - 2 for r in rvals]

plot(rvals, abs.(omvals))
