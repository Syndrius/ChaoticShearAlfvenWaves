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
    
    #this should have a 7/3 island at 0.8 ish
    #and a 8/3 island at 0.85
    a = 1.05
    b = 0.5
    c = 2.35
    return a + b*r^2 + c*r^4, 2 * b * r + 4*c*r^3
    #return a + b*r^2 + c*r^10 + d*r^4, 2*b*r + 10*c*r^9 + 4*d*r^3
    #return a + b*r^6, 6 * b * r^5
    #return a + b*r^4, 4 * b * r^3
end
function dens_prof(r::Float64)
    a = 1.05
    b = 0.4
    c = 5.0
    d = 6

    #this og version doesn't seem to prevent continuum interaction
    #but does isolate the tae's a bit.
    #so it is still nicer
    #rh = 0.7
    #scale = 0.3
    rh = 0.8
    scale = 0.2

    return 1/2 * (1-tanh((r-rh)/scale))

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
rgrid = init_grid(type=:rf, N=150)
θgrid = init_grid(type=:as, start=0, N=10)
#θgrid = init_grid(type=:af, pf=2, N=10)
ζgrid = init_grid(type=:as, start=-3, N=3)
#ζgrid = init_grid(type=:af, pf=-2, N=4)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%

solver = init_solver(nev = 200, target=0.35, prob=prob)

#%%
evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);
#%%

continuum_plot(evals)
continuum_plot(evals, n=-2)

ind  = find_ind(evals, 0.366)

potential_plot(ϕft, grids, ind, label_max=0.2)

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
using Roots
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
#these two pair pretty well for 7/3 and 8/3, 
#however, there still seems to be a significant amount of interaction
#between tae's and continuum, despite frequency gap being open.
function q_prof(r::Float64)
    #we are going to stick with the q-prof
    #n=-3 harmonics for perturbations
    #if we need to look at tae, we will need to change the density profile
    
    #this should have a 7/3 island at 0.8 ish
    #and a 8/3 island at 0.85
    a = 1.05
    b = 0.5
    c = 2.35
    return a + b*r^2 + c*r^4, 2 * b * r + 4*c*r^3
    #return a + b*r^2 + c*r^10 + d*r^4, 2*b*r + 10*c*r^9 + 4*d*r^3
    #return a + b*r^6, 6 * b * r^5
    #return a + b*r^4, 4 * b * r^3
end

function dens_prof(r::Float64)
    a = 1.05
    b = 0.4
    c = 5.0
    d = 6

    rh = 0.7
    scale = 0.3

    return 1/2 * (1-tanh((r-rh)/scale))

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
