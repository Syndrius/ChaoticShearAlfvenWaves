

using MID
using Printf
using DelimitedFiles
using Plots; plotlyjs()

#would be quite nice to have an automatic convergence test function...
#ie just pass in a couple of lists and it does it all??
#considering we will be doing that often?
Nlist = [200, 500, 1000, 1500, 2000, 3000, 4000]
dlist = [-4e-7, -4e-8, -4e-9, -4e-10, -4e-11, -4e-12]
dlabs = [7, 8, 9, 10, 11, 12]

display(split.(string.(dlist), "-"))

dlab = parse(Int64, split(string(dlist[2]), "-")[end])


dir_base = "data/singular_convergence/"

geo = GeoParamsT(R0=20.0)

for N in Nlist
    #rgrid = clustered_grid(N, 0.92, 0.98, 0.2)
    rgrid = clustered_grid(N, 0.75, 0.8, 0.2)
    #grids = init_grids(rgrid=rgrid, mstart=1, mcount=4, nstart=-2, ncount=1);
    grids = init_grids(rgrid=rgrid, mstart=1, mcount=2, nstart=-1, ncount=1);
    for (i, δ) in enumerate(dlist)
        #prob = init_problem(q=Axel_q, geo=geo, δ=δ, dens=axel_dens);
        prob = init_problem(q=singular_bowden_q, geo=geo, δ=δ, dens=bowden_singular_dens);

        dir = dir_base * @sprintf("N_%s_delta_%s/", N, dlabs[i])
        mkdir(dir)

        inputs_to_file(prob=prob, grids=grids, dir=dir)

        ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=(0.33/geo.R0)^2, reconstruct=false);

        eigvals_to_file(ω=ω, filename=dir*"vals.dat")

        eigfuncs_to_file(ϕ=ϕ, filename=dir*"funcs.dat")

    end
end


ωlist = zeros(ComplexF64, length(Nlist), length(dlist))


for (j, N) in enumerate(Nlist)
    #rgrid = clustered_grid(N, 0.92, 0.98, 0.2)
    rgrid = clustered_grid(N, 0.75, 0.8, 0.2)
    #grids = init_grids(rgrid=rgrid, mstart=1, mcount=4, nstart=-2, ncount=1);
    #grids = init_grids(rgrid=rgrid, mstart=1, mcount=2, nstart=-1, ncount=1);
    for (i, δ) in enumerate(dlist)
        #prob = init_problem(q=singular_bowden_q, geo=geo, δ=δ, dens=bowden_singular_dens);

        dir = dir_base * @sprintf("N_%s_delta_%s/", N, dlabs[i])

        vals_file = dir * "vals.dat"
        funcs_file = dir * "funcs.dat"
        
        ω = readdlm(vals_file, ',', ComplexF64)
        #ϕ = reconstruct_phi(readdlm(funcs_file, ',', ComplexF64), length(ω), N, grids.pmd.count, grids.tmd.count)

        #most of these are the first index, some weird stuff is going on later though!
        ωlist[j, i] = ω[1]

        #plot_potential(r=rgrid, ϕ=ϕ, ind=1, pmd=grids.pmd, n=1)

    end
end

#plot_potential(r=rgrid, ϕ=ϕ, ind=2, pmd=grids.pmd, n=1)

#at R0=10
#ours converges to -0.02169
#Em converges to -0.01736
#using percentage difference calculator
#diff of 22.1767%


#at R0=20,
#ours converges to -0.0096555
#Em converges to -0.00896
#diff of 7.47744%
#I think this is enough verification tbh.


p = scatter()#ylimits=(-0.012, -0.008))

for i in 1:length(dlist)
    #scatter!(Nlist, imag.(ωlist[:, i].^2))# ./ real.(ωlist[:, i]))
    scatter!(Nlist, imag.(ωlist[:, i]) ./ real.(ωlist[:, i]))
end
display(p)
display(imag(ωlist[end, end-1]^2))
savefig(p, dir_base*"convergence_zoom.png")









rgrid = clustered_grid(Nlist[end], 0.75, 0.8, 0.2)
grids = init_grids(rgrid=rgrid, mstart=1, mcount=2, nstart=-1, ncount=1);


prob = init_problem(q=singular_bowden_q, geo=geo, δ=dlist[end], dens=bowden_singular_dens);

dir = dir_base * @sprintf("N_%s_delta_%s/", Nlist[end], dlabs[end])

vals_file = dir * "vals.dat"
funcs_file = dir * "funcs.dat"

ω = readdlm(vals_file, ',', ComplexF64)
ϕ = reconstruct_phi(readdlm(funcs_file, ',', ComplexF64), length(ω), Nlist[end], grids.pmd.count, grids.tmd.count)

plot_potential(r=rgrid, ϕ=ϕ, ind=8, pmd=grids.pmd, n=1)

reconstruct_continuum(ϕ=ϕ, ω=ω, grids=grids)