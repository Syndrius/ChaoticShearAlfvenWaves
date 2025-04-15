#testing the jacobain of the qfm transformation.
#in particular, for our chaos base case, using a depth of 6, we get numerical issues. Can be distinguish this from the case with a depth of 5?
using MID
using MIDViz
using JLD2

R0=10.0

#amp needs further thought!
#define the non-resonant island
k = 0.00022
isl = init_island(m0=5, n0=-2, A=k/5)
isl2 = init_island(m0=7, n0=-3, A=k/7)

geo = init_geo(R0=R0)

#to solve non-Hermitian
#flr = MID.Structures.FLRT(δ = 1e-18)
prob = init_problem(q=qfm_benchmark_q, met=:cylinder, geo=geo, isl=isl, isl2=isl2)#, flr=flr)
#so this benchmark case works fine for r<0.4. Unsure why..
k = 0.05
isl = init_island(m0=3, n0=2, A=k/5)
isl2 = init_island(m0=7, n0=-3, A=0.0)
prob = init_problem(q=qfm_benchmark_q, met=:cylinder, geo=geo, isl=isl, isl2=isl2)#, flr=flr)
#this case also has a cooked jacobian, seems to mirror the islands which is interesting...
k = 0.005
isl = init_island(m0=5, n0=-2, A=k/5)
isl2 = init_island(m0=7, n0=-3, A=0.0)
prob = init_problem(q=qfm_benchmark_q, met=:cylinder, geo=geo, isl=isl, isl2=isl2)#, flr=flr)

#%%
qlist, plist = farey_tree(4, 2, 1, 3, 1)
guess_list = 0.5 .* ones(length(qlist));
#needed for the single island chain case, this one is probably at the sepratrix I guess.
deleteat!(qlist, 8)
deleteat!(plist, 8)
deleteat!(qlist, 12)
deleteat!(plist, 12)
@time surfs = construct_surfaces(plist, qlist, guess_list, prob);
plot_surfs(surfs)
#%%
save_object("surfs5.jld2", surfs)
#%%


    #=
    rvals = []
    θvals = []
    ζvals = []

    #Main loop just gets all the actual values that are evaluated.
    for i in 1:grids.r.N-1, j in 1:grids.θ.N, k in 1:grids.ζ.N 

        #takes the local ξ arrays to a global arrays around the grid point.
        r, θ, ζ, dr, dθ, dζ = MID.Basis.local_to_global(i, j, k, ξr, ξθ, ξζ, rgrid, θgrid, ζgrid) 

        for x1 in r, x2 in θ, x3 in ζ
            push!(rvals, x1)
            push!(θvals, x2)
            push!(ζvals, x3)

        end

    end
    =#
#%%
function compute_jac(prob, grids, surfs)

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
    ξr, wgr = MID.Construct.FastGaussQuadrature.gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = MID.Construct.FastGaussQuadrature.gausslegendre(grids.θ.gp)
    ξζ, wgζ = MID.Construct.FastGaussQuadrature.gausslegendre(grids.ζ.gp)

    #struct for storing the intermediate data for the coordinate transform
    CT = MID.QFM.CoordTsfmT()
    #jac = zeros(length(rvals), length(θvals), length(ζvals))
    jac = zeros(grids.r.N, grids.θ.N, grids.ζ.N)
    jac_tor = zeros(grids.r.N, grids.θ.N, grids.ζ.N)

    #for (i, r) in enumerate(rvals), (j, θ) in enumerate(θvals), (k, ζ) in enumerate(ζvals)
    for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)
        MID.QFM.coord_transform!(r, θ, ζ, CT, surf_itp, sd)
        MID.Geometry.toroidal_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], 4.0)
        MID.Equilibrium.compute_B!(tor_B, tor_met, prob.q, prob.isl, prob.isl2, CT.coords[1], CT.coords[2], CT.coords[3])
        MID.QFM.met_transform!(tor_met, qfm_met, CT)
        MID.QFM.B_transform!(tor_B, qfm_B, qfm_met, CT)

        jac[i, j, k] = qfm_met.J[1]
        jac_tor[i, j, k] = tor_met.J[1]
        #jac[i, j, k] = qfm_B.B[1]
        #jac_tor[i, j, k] = qfm_B.B[2]
    end
    return jac, jac_tor
end
#%%
#this grid caused issues.
Nr = 100
Nθ = 20
Nζ = 1
#ok so the jacobian below ~r=0.4 is completly cooked. No idea why. unsure how to fix.
#I guess this is because of the axis? not really sure how to fix that though?
rgrid = init_grid(type=:rf, N = Nr, start = 0.2, stop =0.7)
θgrid = init_grid(type=:af, N = Nθ, pf=5)
ζgrid = init_grid(type=:af, N = Nζ, pf=-2)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%
surfs5 = load_object("surfs5.jld2");
surfs6 = load_object("surfs6.jld2");
#%%
jac, jac_tor = compute_jac(prob, grids, surfs);
#%%
jac5 = compute_jac(prob, grids, surfs5);
jac6 = compute_jac(prob, grids, surfs6);
#%%
#jacobian for r<0.4 seems to be breaking down, correlated with original jacobian being less than 1.
#so the jacobian is clearly just very large inside the islands, We can perhaps avoid this problemo by just considering the region away from the islands. At least no jacobian value is negative or anything stupid, just very large.
#perhaps we could remove some of the surfaces around the islands.
rgrid_plot, θgrid_plot, _ = MID.Structures.inst_grids(grids);
contourf(θgrid_plot, rgrid_plot, jac[:, :, 1], levels=50)
contourf(θgrid_plot, rgrid_plot, jac_tor[:, :, 1], levels=50)
contourf(θgrid_plot, rgrid_plot, jac_tor[:, :, 1] ./ jac[:, :, 1], levels=50)
contourf(θgrid_plot, rgrid_plot, jac5[:, :, 3])
contourf(θgrid_plot, rgrid_plot, jac6[:, :, 3])

