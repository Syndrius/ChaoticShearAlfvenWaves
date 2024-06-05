using MID
#we will probably remove the other damping benchmarks eventually.

#This benchmark compares to Bowden and Hole 2015.
#Here we compute the continuum damping of a tae.
#This paper solves the simplified two mode tae equation derived by Berk et al 1992
#This equation makes many approximation about large aspect ratio's that we do not.
#We have implemented our own simple finite element solver directly for this equation to compare to.
#we start with R=10, as done in the paper.

#may need to make sure these are the actual converged values.
N = 3000
δ = -4e-9
R0 = 10.0

tae_freq = (0.3259/R0)^2


grids = init_grids(N=N, sep1=0.75, sep2=0.8, frac=0.2, mstart=1, mcount=2, nstart=-1, ncount=1);

geo = GeoParamsT(R0=R0)

prob = init_problem(q=singular_bowden_q, geo=geo, δ=δ, dens=bowden_singular_dens);

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq, reconstruct=true);


rgrid = MID.Misc.clustered_grid(N, 0.75, 0.8, 0.2);
#this function computes the damping directly from Berk's equation.
ω_two, Em, Em1 = two_mode(rgrid=rgrid, R0=R0, m=1, n=-1, δ=δ, σ=tae_freq);

display(ω_two[1])
#our code predictions a damping ratio of ~-0.022
dr = imag(ω[1])/real(ω[1])
#Berks equation gives a damping ratio of ~-0.0174, giving almost identical results to the quoted 
#-0.0175 in Bowden and Hole.
dr_two = imag(ω_two[1])/real(ω_two[1])

#our result is ~22% different
display(abs(dr-dr_two)/((dr + dr_two)/2) * 100)

#however, if we increase the aspect ratio so that the equation are more similar,
R0 = 20.0
geo = GeoParamsT(R0=R0)
tae_freq = (0.332/R0)^2 #tae frequency shifts a bit with the new aspect ratio

prob = init_problem(q=singular_bowden_q, geo=geo, δ=δ, dens=bowden_singular_dens);

ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=tae_freq, reconstruct=false);
ω_two, Em, Em1 = two_mode(rgrid=rgrid, R0=R0, m=1, n=-1, δ=δ, σ=tae_freq);

display(ω_two[1])
#our damping ratio is now ~-0.000966
dr = imag(ω[1])/real(ω[1])
#And the damping ration from Berk's equation is ~-0.00896
dr_two = imag(ω_two[1])/real(ω_two[1])

#yeilding a potential difference of only ~7.5%
#and we would expect in the limit R0→∞ or ϵ→0, these two results would converge.
display(abs(dr-dr_two)/((dr + dr_two)/2) * 100)

#######################################################

"""

This section is for running the convergence tests to verify the balues computed above.

"""

#=

using Plots
using LaTeXStrings

#alternativly, this can be seen with a proper convergence test,
#start with R0=10 case.
Nlist = [200, 500, 1000, 1500, 2000, 3000, 4000]
#probably shouldn't include δ=-4e-12 tbh
δlist = [-4e-7, -4e-8, -4e-9, -4e-10, -4e-11, -4e-12]
R0 = 10.0

tae_freq = (0.3259/R0)^2
grids = init_grids(N=10, sep1=0.75, sep2=0.8, frac=0.2, mstart=1, mcount=2, nstart=-1, ncount=1);
geo = GeoParamsT(R0=R0)

prob = init_problem(q=singular_bowden_q, geo=geo, δ=1.0, dens=bowden_singular_dens);
convergence_test(grids, prob, Nlist, δlist, "data/damping_convergence/", tae_freq)
two_mode_convergence(Nlist, δlist, "data/two_mode_damping/", R0, tae_freq)

#then we consider the R0=20 case.
R0 = 20.0
tae_freq = (0.332/R0)^2
grids = init_grids(N=10, sep1=0.75, sep2=0.8, frac=0.2, mstart=1, mcount=2, nstart=-1, ncount=1);
geo = GeoParamsT(R0=R0)

prob = init_problem(q=singular_bowden_q, geo=geo, δ=1.0, dens=bowden_singular_dens);
convergence_test(grids, prob, Nlist, δlist, "data/damping_convergence_20/", tae_freq)
two_mode_convergence(Nlist, δlist, "data/two_mode_damping_20/", R0, tae_freq)





ωlist = read_convergence_data(Nlist, δlist, "data/damping_convergence/");
ωlist_2m = read_convergence_data(Nlist, δlist, "data/two_mode_damping/");

ωlist_20 = read_convergence_data(Nlist, δlist, "data/damping_convergence_20/");
ωlist_2m_20 = read_convergence_data(Nlist, δlist, "data/two_mode_damping_20/");




shapes = [:circle, :square, :dtriangle, :pentagon, :star, :diamond]
p = scatter(xlabel="Grid Points", ylabel=L"Im($\Omega$)/Re($\Omega$)", left_margin=4Plots.mm, yguidefontsize=16, xguidefontsize=16, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=11, ylimits=(-0.023, -0.020))
for i in 1:1:length(δlist)
    δlab = δlist[i]
    scatter!(Nlist, imag.(ωlist[:, i]) ./ real.(ωlist[:, i]), markershape=shapes[i], label=L"\delta=%$δlab", opacity=0.8, markersize=6)
end
conv_val = imag.(ωlist[6, 4]) ./ real.(ωlist[6, 4])
hline!(conv_val .* ones(length(Nlist)), linestyle=:dash, color=:black, label=L"\Omega=%$(round(conv_val; digits=4))")
display(p)
#savefig(p, "data/damping_convergence.png")
#savefig(p, "data/damping_convergence_zoom.png")


p = scatter(xlabel="Grid Points", ylabel=L"Im($\Omega$)/Re($\Omega$)", left_margin=4Plots.mm, yguidefontsize=16, xguidefontsize=16, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=11)#, ylimits=(-0.0185, -0.0165))
for i in 1:1:length(δlist)
    δlab = δlist[i]
    scatter!(Nlist, imag.(ωlist_2m[:, i]) ./ real.(ωlist_2m[:, i]), markershape=shapes[i], label=L"\delta=%$δlab", opacity=0.8)
end
conv_val = imag.(ωlist_2m[6, 4]) ./ real.(ωlist_2m[6, 4])
hline!(conv_val .* ones(length(Nlist)), linestyle=:dash, color=:black, label=L"\Omega=%$(round(conv_val; digits=4))")
display(p)
#savefig(p, "data/two_mode_damping.png")
#savefig(p, "data/two_mode_damping_zoom.png")


p = scatter(xlabel=LaTeXString("Grid Points"), ylabel=L"Im($\Omega$)/Re($\Omega$)", left_margin=4Plots.mm, yguidefontsize=16, xguidefontsize=16, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)#, ylimits=(-0.0105, -0.009))
for i in 1:1:length(δlist)
    δlab = δlist[i]
    scatter!(Nlist, imag.(ωlist_20[:, i]) ./ real.(ωlist_20[:, i]), markershape=shapes[i], label=L"\delta=%$δlab", opacity=0.8)
end
conv_val = imag.(ωlist_20[6, 4]) ./ real.(ωlist_20[6, 4])
hline!(conv_val .* ones(length(Nlist)), linestyle=:dash, color=:black, label=L"\Omega=%$(round(conv_val; digits=4))")
display(p)
#savefig(p, "data/damping_convergence_20.png")
#savefig(p, "data/damping_convergence_20_zoom.png")


p = scatter(xlabel="Grid Points", ylabel=L"Im($\Omega$)/Re($\Omega$)", left_margin=4Plots.mm, yguidefontsize=16, xguidefontsize=16, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=11)#, ylimits=(-0.0095, -0.008))
for i in 1:1:length(δlist)
    δlab = δlist[i]
    scatter!(Nlist, imag.(ωlist_2m_20[:, i]) ./ real.(ωlist_2m_20[:, i]), markershape=shapes[i], label=L"\delta=%$δlab", opacity=0.8)
end
conv_val = imag.(ωlist_2m_20[6, 4]) ./ real.(ωlist_2m_20[6, 4])
hline!(conv_val .* ones(length(Nlist)), linestyle=:dash, color=:black, label=L"\Omega=%$(round(conv_val; digits=4))")
display(p)
#savefig(p, "data/two_mode_damping_20.png")
#savefig(p, "data/two_mode_damping_20_zoom.png")

=#