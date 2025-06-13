
#why are our matrices not perfectly Hermitian?

using MID
using LinearAlgebra
using JLD2
#%%

rgrid = init_grid(type=:rf, N = 3, gp=6)
θgrid = init_grid(type=:af, N = 2, gp=6, pf=2)
ζgrid = init_grid(type=:af, N = 2, gp=6, pf=-1)
#θgrid = init_grid(type=:as, N = 2, start=1)
#ζgrid = init_grid(type=:as, N = 2, start=-2)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%
surfs = load_object("/Users/matt/phd/MID/data/surfaces/island_21a.jld2");
#%%
geo = init_geo(R0=10.0)

isl21a = init_island(m0=2, n0=-1, w=0.05, r0=0.5, qp=2.0)
#start with no islands
prob = init_problem(geo=geo, q=MID.Equilibrium.island_equiv_q, met=:cylinder, isl=isl21a)

#%%

W, I, arr = MID.Construct.construct(prob, grids)
W, I = MID.Construct.construct(prob, grids)
#%%
Wmat = Matrix(W);
Imat = Matrix(I);
W2mat = Matrix(W);
I2mat = Matrix(I);

ishermitian(Wmat)
ishermitian(Imat)
display(Wmat[5:10, 5:10])
display(Imat[5:10, 5:10])
#%%
#this is less symmetric lol
#guess the way they combine the data is inconsistent?
#Petsc has the option to just ignore the lower triangle completly, which may be our best bet!
#i am assuming the sparse matrix is constructed by iterating through all the rows first, 
#so floating point error is expected.
display(maximum(abs.(Wmat .- Wmat')))
display(maximum(abs.(Imat .- Imat')))
#this should be zero if the only issue is the order that integration is done.
#tis close for I, but not exact
display(Imat .- I2mat')
#%%
isapprox(Imat, Imat', rtol=1e-17)
isapprox(Wmat, Wmat', rtol=1e-19)

display(arr)
display(arr2)
sum(arr)
sum(arr2)

#ok this means that each individual value is not the same...
#so problem is not occuring with the sum here.
sort(real.(arr)) .- sort(real.(arr2))
sum(arr2)
sort(real.(arr))
sort(real.(arr2))

#ok so looks like there are multiple problemos
#W from Tl is actually not symmetric, this may be quite a serious issue, however, this could just be due to rounding
#Each Gauss integrate is not giving exactly the same, even for I. 
#this means we cannot explain the symmetry loss by addition of floats. This seems to be a very small issue


C = [[0.0  -0.0551168   -7.67239e-5  0.0  0.0554179  0.0997434  0.0        0.0        0.0] ;[ 0.0  -2.74478e-6  -4.94017e-6  0.0  0.0        0.0        0.0554179  0.0997434  0.0] ; [0.0  -0.0         -0.0         0.0  0.0        0.0        0.0        0.0554179  0.0997434]]
D = [[1.04272    0.081085    -0.0]; [0.081085   1.05548     -0.00485749]; [0.0       -0.00485749   2.24893e-5]]


res = (C' * D) * C
res .- res'
ishermitian((C' * D) * C)

res2 = zeros(9, 9)

for i in 1:9, j in 1:9, k in 1:3, l in 1:3
    res2[i, j] = C[k, i] * D[k, l] * C[l, j]
end

ishermitian(res2)
