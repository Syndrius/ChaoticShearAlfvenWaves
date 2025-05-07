
#trying to get our chaotic region to be better, by adding in more than 2 islands.

using MID
using MIDViz
using JLD2
using Plots; plotlyjs()
using Plots; gr()
#%%
function rationals(depth, min, max)

    rats = Tuple[]
    for i in 1:depth
        j = 1
        while (j / i < max)
            if j/i > min
                a = gcd(j, i)
                if !((j, i) in rats)
                    if a==1
                        push!(rats, (j, i))
                    elseif !((Int(j/a), Int(i/a)) in rats)
                        push!(rats, (Int(j/a), Int(i/a)))
                    end
                end
            end
            j += 1
        end
    end
    return rats
end


function q_prof(r::Float64)
    a = 4/3
    b = 25/24
    #a = 359/360
    #b = 10/9
    a = 223/195
    b = 40/39
    a = 254/255
    b = 80/51
    a = 21/20
    b = 1.0
    a = 1.05
    b = 0.45
    return a + b * r^2, 2 * b * r
end
#%%
q_prof(0.05)
rats = rationals(7, q_prof(0.0)[1], q_prof(1.0)[1])
rats = rationals(15, 1.1, 1.23)
vals = sort([i[1]/i[2] for i in rats])
[i[1]/i[2] for i in rats]
#%%
geo = init_geo(R0=4.0)

k = 0.0003
isl1 = init_island(m0=7, n0=-5, A=0.0*k/12)
isl2 = init_island(m0=10, n0=-7, A=0.0*k/10)
isl3 = init_island(m0=5, n0=-3, A=1.0*k/7)
isl4 = init_island(m0=7, n0=-4, A=0.0*k/7)
isl5 = init_island(m0=8, n0=-5, A=1.0*k/3)
isl6 = init_island(m0=9, n0=-5, A=1.0*k/9)
isl7 = init_island(m0=11, n0=-7, A=0.0*k/7)
isl8 = init_island(m0=13, n0=-9, A=0.0*k/8)
#I think we want to do whatever is needed to not have any m or n < ~5
#these cause way to much deformation of the flux surfaces.
isl9 = init_island(m0=3, n0=-2, A=1.0*k/12)

#isls = [isl1, isl2, isl3, isl4, isl5, isl6, isl7, isl8, isl9]
isl1 = init_island(m0=13, n0=-8, A=1.0*k/13)
isl2 = init_island(m0=9, n0=-5, A=1.0*k/9)
isl3 = init_island(m0=12, n0=-7, A=1.0*k/12)
isl4 = init_island(m0=5, n0=-3, A=1.0*k/5)
isl5 = init_island(m0=8, n0=-5, A=0.0*k/11)
isl6 = init_island(m0=7, n0=-4, A=1.0*k/7)
isl7 = init_island(m0=11, n0=-7, A=1.0*k/11)

isls = [isl1, isl2, isl3, isl4, isl5, isl6, isl7]


isl1 = init_island(m0=6, n0=-5, A=1.0*k/6)
isl2 = init_island(m0=11, n0=-9, A=1.0*k/11)
isl3 = init_island(m0=9, n0=-7, A=1.0*k/6)
isl4 = init_island(m0=17, n0=-13, A=1.0*k/17)
isl5 = init_island(m0=13, n0=-11, A=1.0*k/13)
isl6 = init_island(m0=5, n0=-4, A=1.0*k/5)
isl7 = init_island(m0=13, n0=-10, A=1.0*k/13)
isl8 = init_island(m0=14, n0=-11, A=1.0*k/14)

isls = [isl1, isl2, isl3, isl4, isl5, isl6, isl7, isl8]

prob = init_problem(geo=geo, q=q_prof, isls=isls)
unprob = init_problem(geo=geo, q=q_prof)
#%%

rgrid = init_grid(type=:rf, start=0.05, stop=0.95, N=100)
θgrid = init_grid(type=:as, start=0, N=6)
ζgrid = init_grid(type=:as, start=-4, N=4)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%
solver = init_solver(targets=[0.2, 0.3, 0.4], nev=200, prob=unprob)
#%%
evals_norm, ϕ_norm, ϕft_norm = compute_spectrum(prob=unprob, grids=grids, solver=solver);

continuum_plot(evals_norm)

ind_norm = find_ind(evals_norm, 0.342)

potential_plot(ϕft_norm, grids, ind_norm, label_max=0.5)
#%%
Ntraj = 100;
#rlist = collect(LinRange(0.4, 0.8, Ntraj));
rlist = collect(LinRange(0.1, 1.0, Ntraj));
Nlaps = 500;

poincare_plot(prob, Nlaps, Ntraj, rlist)
#%%
#ok lets get some god damn surfaces then.
q_prof(0.0)
q_prof(1.0)
rats1 = rationals(7, q_prof(0.0)[1], q_prof(1.0)[1])
ql1 = [i[1] for i in rats1]
pl1 = [i[2] for i in rats1]
gl1 = [sqrt(i[1]/i[2] - q_prof(0.0)[1])/sqrt(q_prof(1.0)[1] - q_prof(0.0)[1]) for i in rats1]
#increasing the numbers did have an effect, surfaces no longer overlap
#there are still jagged peices of garbage though
@time surfs1 = construct_surfaces(pl1, ql1, gl1, prob, MM=4, M=32, N=12);
plot_surfs(surfs1)
length(rats1)
length(surfs1)
#%%
rats2 = [(11, 10), (10, 9), (9, 8), (11, 8), (11, 9)]
ql2 = [i[1] for i in rats2]
pl2 = [i[2] for i in rats2]
gl2 = [sqrt(i[1]/i[2] - q_prof(0.0)[1])/sqrt(q_prof(1.0)[1] - q_prof(0.0)[1]) for i in rats2]
@time surfs2 = construct_surfaces(pl2, ql2, gl2, prob, MM=4, M=32, N=12);
plot_surfs(surfs2)
#%%
rats3 = [(29, 20), (13, 10), (15, 11), (20, 19)]
ql3 = [i[1] for i in rats3]
pl3 = [i[2] for i in rats3]
gl3 = [sqrt(i[1]/i[2] - q_prof(0.0)[1])/sqrt(q_prof(1.0)[1] - q_prof(0.0)[1]) for i in rats3]
@time surfs3 = construct_surfaces(pl3, ql3, gl3, prob, MM=4, M=32, N=12);
plot_surfs(surfs3)
#%%
rats4 = [(15, 14), (17, 16), (18, 17)]
ql4 = [i[1] for i in rats4]
pl4 = [i[2] for i in rats4]
gl4 = [sqrt(i[1]/i[2] - q_prof(0.0)[1])/sqrt(q_prof(1.0)[1] - q_prof(0.0)[1]) for i in rats4]
@time surfs4 = construct_surfaces(pl4, ql4, gl4, prob, MM=4, M=32, N=12);
plot_surfs(surfs4)
#%%
rats5 = [(13, 12)]
ql5 = [i[1] for i in rats5]
pl5 = [i[2] for i in rats5]
gl5 = [sqrt(i[1]/i[2] - q_prof(0.0)[1])/sqrt(q_prof(1.0)[1] - q_prof(0.0)[1]) for i in rats5]
@time surfs5 = construct_surfaces(pl5, ql5, gl5, prob, MM=4, M=32, N=12);
plot_surfs(surfs5)
#%%
rats6 = [(14, 11), (19, 18), (12, 11)]
ql6 = [i[1] for i in rats6]
pl6 = [i[2] for i in rats6]
gl6 = [sqrt(i[1]/i[2] - q_prof(0.0)[1])/sqrt(q_prof(1.0)[1] - q_prof(0.0)[1]) for i in rats6]
@time surfs6 = construct_surfaces(pl6, ql6, gl6, prob, MM=4, M=32, N=12);
plot_surfs(surfs6)
#%%
rats7 = [(21, 19), (13, 11), (15, 13), (17, 13)]
ql7 = [i[1] for i in rats7]
pl7 = [i[2] for i in rats7]
gl7 = [sqrt(i[1]/i[2] - q_prof(0.0)[1])/sqrt(q_prof(1.0)[1] - q_prof(0.0)[1]) for i in rats7]
@time surfs7 = construct_surfaces(pl7, ql7, gl7, prob, MM=4, M=32, N=12);
plot_surfs(surfs7)
#%%
curr_surfs = surfs1;
#not sure 6 actually helps.
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5, surfs6);
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5);
curr_surfs = vcat(surfs1, surfs2, surfs3, surfs4, surfs5, surfs6, surfs7);

plot_surfs(curr_surfs)
#(9, 5) is ind 9.
#12, 7 is 25
#11, 6, is 19
plot_surfs(surfs1[9:9])
display(rats1)
display(length(rats1))
curr_surfs = vcat(surfs1[1:13], surfs1[15:18], surfs1[20:24], surfs1[26:end]);
plot_surfs(curr_surfs)
save_object("total_low_shear_surfs.jld2", curr_surfs)
curr_surfs = load_object("total_low_shear_surfs.jld2");
#%%

evals, ϕ, ϕft = compute_spectrum_qfm(grids=grids, prob=prob, solver=solver, surfs=curr_surfs);

ontinuum_plot(evals)

ind = find_ind(evals, 0.343)

potential_plot(ϕft, grids, ind, label_max=0.5)
#%%
function compute_flux(prob, grids, surfs)

    #instantiate the grids into arrays. 
    rgrid, θgrid, ζgrid = MID.Structures.inst_grids(grids)

    #initialise the two structs to store the metric and the magnetic field.
    tor_met = MID.Geometry.MetT()
    qfm_met = MID.Geometry.MetT()
    tor_B = MID.Equilibrium.BFieldT()
    qfm_B = MID.Equilibrium.BFieldT()

    #creates the interpolations for the surfaces.
    surf_itp, sd = MID.QFM.create_surf_itp(surfs)

    #compute the gaussian qudrature points for finite elements.
    #ξr, wgr = MID.Construct.FastGaussQuadrature.gausslegendre(grids.r.gp) #same as python!
    #ξθ, wgθ = MID.Construct.FastGaussQuadrature.gausslegendre(grids.θ.gp)
    #ξζ, wgζ = MID.Construct.FastGaussQuadrature.gausslegendre(grids.ζ.gp)

    #struct for storing the intermediate data for the coordinate transform
    CT = MID.QFM.CoordTsfmT()
    #jac = zeros(length(rvals), length(θvals), length(ζvals))
    jac = zeros(length(rgrid), length(θgrid), length(ζgrid))
    djac = zeros(3, length(rgrid), length(θgrid), length(ζgrid))
    B = zeros(3, length(rgrid), length(θgrid), length(ζgrid))

    #for (i, r) in enumerate(rvals), (j, θ) in enumerate(θvals), (k, ζ) in enumerate(ζvals)
    for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)
        MID.QFM.coord_transform!(r, θ, ζ, CT, surf_itp, sd)
        MID.Geometry.toroidal_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], 4.0)
        MID.Equilibrium.compute_B!(tor_B, tor_met, prob.q, prob.isls, CT.coords[1], CT.coords[2], CT.coords[3])
        MID.QFM.met_transform!(tor_met, qfm_met, CT)
        MID.QFM.B_transform!(tor_B, qfm_B, qfm_met, CT)

        jac[i, j, k] = qfm_met.J[1]
        djac[:, i, j, k] = qfm_met.dJ[:]
        B[:, i, j, k] = qfm_B.B[:]
        #jac_tor[i, j, k] = tor_met.J[1]
        #jac[i, j, k] = qfm_B.B[1]
        #jac_tor[i, j, k] = qfm_B.B[2]
    end
    return B, jac, djac
end
#%%
B, jac, djac = compute_flux(prob, grids, curr_surfs);
#%%
rgrid_plot, θgrid_plot, _ = MID.Structures.inst_grids(grids);
contourf(θgrid_plot, rgrid_plot, jac[:, :, 1], levels=50, title="Jacobian")
contourf(θgrid_plot, rgrid_plot, djac[1, :, :, 1], levels=50, title="dJdr")
contourf(θgrid_plot, rgrid_plot, djac[2, :, :, 1], levels=50, title="dJdt")
#contourf(θgrid_plot, rgrid_plot, B[1, :, :, 1], levels=50, title="B^s")
contourf(θgrid_plot, rgrid_plot, (B[1, :, :, 1]).^2, levels=50, title="B^s^2")
#%%
B2mean = zeros(length(surfs_list));
djacdrmean = zeros(length(surfs_list));
djacdθmean = zeros(length(surfs_list));

for i in 1:length(surfs_list)
    curr_surfs = surfs_list[i]

    B, jac, djac = compute_flux(prob, grids, curr_surfs);

    B2mean[i] = mean(B[1, :, :, :] .^ 2)
    djacdrmean[i] = mean(abs.(djac[1, :, :, :]))
    djacdθmean[i] = mean(abs.(djac[2, :, :, :]))
end
#%%
