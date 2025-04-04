
using MID
using Plots#; plotlyjs()
using LaTeXStrings
using Printf

#want to create an overlay plot like Axel's case, where we use the continuum function to find the continuum, but overlay the local area around the tae.


N = 100;
grids_cont = init_grids(N=N, mstart=2, mcount=2, nstart=-2, ncount=1);

geo = GeoParamsT(R0=10.0)
prob = init_problem(q=island_damping_q, geo=geo); 


ω_cont = continuum(prob=prob, grids=grids_cont);
#plot_continuum(ω = ω_cont, grids=grids)


#ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=true);


#reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids, ymax=0.75)

#maybe should start making the target freq the normalised version and automatically un-normalise it?
#would prevent doing this shite all the time.
tae_freq = (0.382807^2/geo.R0^2)



N = 2000;
grids = init_grids(N=N, mstart=2, mcount=2, nstart=-2, ncount=1);


ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq, nev=1);


#reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)
plot_potential(grids=grids, ϕ=ϕ, ind=1, n=1)

#cureently best sol seems to be having nev=1, to just get the tae.
overlay_plot(ω_cont, grids_cont, ω, ϕ, grids)#, "island_damping_overlay.png")


#basic functionality works, bit clunky, bit garbage tbh.
function overlay_plot(ω_cont, grids_cont, ω, ϕ, grids, filename=nothing)

    p = scatter(xlabel=L"r", ylabel=L"\frac{\omega  R_0}{v_A}", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)#, guidefontvalign=:hcentre)
    mlist_cont = (grids_cont.pmd.start:grids_cont.pmd.incr:grids_cont.pmd.start + grids_cont.pmd.incr * grids_cont.pmd.count)[1:end-1]
    rgrid = collect(LinRange(0, 1, grids_cont.rd.N))[2:end]
    n = -2
    for (i, m) in enumerate(mlist_cont)
        scatter!(rgrid, ω_cont[:, i], label=@sprintf("(%d, %d)", m, n))#, linestyle=:dash, linewidth=4)
    end


    omdata = zeros(length(ω)) 
    rdata = zeros(length(ω))
    #col = zeros(length(ω))

    col = Tuple{Int, Int}[]

    #rm = zeros(Int, pmd.count, tmd.count)
    #ϕm = zeros(pmd.count, tmd.count)
    labels = zeros(grids.pmd.count, grids.tmd.count)

    #may be better to use other modules and call spectral_grid or whatever
    mlist = (grids.pmd.start:grids.pmd.incr:grids.pmd.start + grids.pmd.incr * grids.pmd.count)[1:end-1]
    nlist = (grids.tmd.start:grids.tmd.incr:grids.tmd.start + grids.tmd.incr * grids.tmd.count)[1:end-1]
    
    rgrid = MID.construct_rgrid(grids)
    #labels = [collect(pmd.start:pmd.start+pmd.count), collect(tmd.start:tmd.start+tmd.count)]

    for i in 1:1:length(ω)
        #still assumes only a single n value!
        rm = zeros(Int, grids.pmd.count, grids.tmd.count)
        ϕm = zeros(grids.pmd.count, grids.tmd.count)
        for j in 1:grids.pmd.count

            for k in 1:grids.tmd.count
                rm[j, k] = argmax(abs.(real.(ϕ[i, :, j, k])))
                ϕm[j, k] = abs.(real.(ϕ[i, rm[j, k], j, k]))
            end
        end

        max_mode = argmax(ϕm)

        #display(max_mode)
        #display(max_mode[2])
        #display(rm[max_mode])

        rdata[i] = rgrid[rm[max_mode]]
        #col[i] = (mlist[max_mode[1]], nlist[max_mode[2]])
        push!(col, (mlist[max_mode[1]], nlist[max_mode[2]]))
        omdata[i] = abs.(ω[i]) #already normalised now!
        #omdata[i] = sqrt.(abs.(ω[i])) * R0


    end

    scatter!(rdata, omdata, group=col, label=false)
    display(p)

    if !isnothing(filename)
        savefig(p, filename)
    end


end
