
#can we get extrapolation of qfm surfaces to work.
#first we need to understand exactly what we are interpolating.
using MID
using MIDViz
using Plots
#%%

R0=4.0

#amp needs further thought!
#define the non-resonant island
k = 0.005
isl = init_island(m0=3, n0=2, A=k)
#isl2 = init_island(m0=7, n0=-3, A=k/7)

geo = init_geo(R0=R0)

#to solve non-Hermitian
#flr = MID.Structures.FLRT(δ = 1e-18)
prob = init_problem(q=qfm_benchmark_q, geo=geo, isl=isl)#, isl2=isl2)#, flr=flr)
#%%
qlist, plist = farey_tree(4, 2, 1, 3, 1)


guess_list = 0.5 .* ones(length(qlist));
guess_list[1] = 0.2
guess_list[2] = 0.8

#%%
#For depth of 3, giving 9 surfaces, took ~70s
#fkn stupid af this is plist then qlist, opposite to farey_tree.
#depth of 5 took 33 mins lol. But basically only 2 surfs where the problemo.
#and they are cooked implying a proper solution was not found
@time surfs = construct_surfaces(plist, qlist, guess_list, prob);
#%%

surf_itp = MID.QFM.create_surf_itp(surfs);

#%%
function get_r_θ(surf_itp, s, ϑ, φ)

    pqMpol = surf_itp.M
    pqNtor =surf_itp.N
    dim1 = pqMpol+1
    dim2 = 2*pqNtor + 1

    mlist = collect(range(0, pqMpol));

    #collect(-pqNtor:0)
    nlist = [collect(0:pqNtor);collect(-pqNtor:-1)] ;

    α = zeros((length(mlist), length(nlist)));

    for i in 1:1:length(mlist)
        for j in 1:1:length(nlist)
            α[i, j] = mlist[i] * ϑ - nlist[j] * φ
        end
    end

    cosα = cos.(α);
    sinα = sin.(α);

    rcos, θsin = MID.QFM.itp_mat(surf_itp, s, deriv=0);
    drcosds, dθsinds = MID.QFM.itp_mat(surf_itp, s, deriv=1)
    display(drcosds)
    display(dθsinds)
    r = 0.0
    θ = ϑ

    for i in 1:dim1, j in 1:dim2
        r += rcos[i, j] * cosα[i, j]
        θ += θsin[i, j] * sinα[i, j]
    end
    return r, θ
end

#%%
s = 0.1
ϑ = 0.2
φ = 0.0

r, θ = get_r_θ(surf_itp, s, ϑ, φ)
#%%

sgrid = LinRange(0.4, 0.6, 10)
sgrid = [0.05]
ϑgrid = LinRange(0, 2π, 50)
rres = zeros(length(sgrid), length(ϑgrid));
θres = zeros(length(sgrid), length(ϑgrid));
for (i, sval) in enumerate(sgrid), (j, ϑval) in enumerate(ϑgrid)

    r, θ = get_r_θ(surf_itp, sval, ϑval, φ)
    rres[i, j] = r
    θres[i, j] = θ
end
#%%
    
#seems like exterpolation should be working. Unclear where the linear algebra exception was coming from!
plot(rres[1, :])
plot(θres[1, :])


#%%
CT = MID.QFM.CoordTsfmT()
MID.QFM.coord_transform!(0.1999, ϑ, φ, CT, surf_itp)
display(CT.JM)




#%%

#see if we can get a minimum working example of derivative interpolations not working
#so the extrapolation doesn't work for derivatives. This explains the problemo, we may need to swap packages. Although ideally we won't.
using BSplineKit

xs = 0.2:0.2:1.2
ys = 2 * cospi.(xs)
S = interpolate(xs, ys, BSplineOrder(6))
ext_S = extrapolate(S, Smooth())
dS = Derivative(1)*S
ext_dS = extrapolate(dS, Smooth())

#ok so this works!
#doesn't seem like it will be very accurate for long distance extrapolations, although maybe this
#is a symptom of this scenario.
#so I think we will still need to be cautious of the extrapolation, but it least it won't error straight away. but a long way away from the surfaces will cause problemos though.
(Derivative(1)*ext_S)(1.3)
ext_dS(1.3)

