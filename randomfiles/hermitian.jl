
#why are our matrices not perfectly Hermitian?

using MID
using LinearAlgebra
using JLD2
#%%

rgrid = init_grid(type=:rf, N = 4, gp=4)
θgrid = init_grid(type=:af, N = 2, gp=6, pf=2)
ζgrid = init_grid(type=:af, N = 2, gp=6, pf=-1)
θgrid = init_grid(type=:as, N = 2, start=1)
ζgrid = init_grid(type=:as, N = 1, start=-1)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%
surfs = load_object("/Users/matt/phd/MID/data/surfaces/island_21a.jld2");
#%%
geo = init_geo(R0=4.0)

isl21a = init_island(m0=2, n0=-1, w=0.05, r0=0.5, qp=2.0)
#isl21a = init_island(m0=2, n0=-1, A=0.001)
#start with no islands
prob = init_problem(geo=geo, q=MID.Equilibrium.island_equiv_q, met=:cylinder, isl=isl21a)

prob = init_problem(q=fu_dam_q, geo=geo)

#%%

#W, I, arr = MID.Construct.construct(prob, grids)
W, I = MID.Construct.construct(prob, grids)
W, I = MID.Construct.construct(prob, grids, surfs)
#%%
Wmat = Matrix(W);
Imat = Matrix(I);
W2mat = Matrix(W);
I2mat = Matrix(I);

eigen(real.(Imat))
eigen(Imat)
isposdef(real.(Imat))

symI = ( Imat .+ Imat') ./ 2

isposdef(symI)
#exact same eigen spectrum as before
#this tells us that our matrices are essentially positive definite.
#need to figure out cases where slepc thinks it is not.
#unsure what to do then
#same seems to be true for qfm cases, although it is defs less posdefinite.
eigen(symI)

#so why is this allowed if there are not positive def?
#and not Hermitian? There must be a threshold, ideally we should set that in Slepc.
eigen(Hermitian(Wmat), Hermitian(Imat))

ishermitian(Wmat)
ishermitian(Imat)
display(Wmat[10:12, 8:10])
display(Imat[6:10, 6:10])
println(Wmat[7, 11])
println(Wmat[11, 7])
#%%
#this is less symmetric lol
#guess the way they combine the data is inconsistent?
#Petsc has the option to just ignore the lower triangle completly, which may be our best bet!
#i am assuming the sparse matrix is constructed by iterating through all the rows first, 
#so floating point error is expected.
display(maximum(abs.(Wmat .- Wmat')))
display(argmax(abs.(Wmat .- Wmat')))
display(maximum(abs.(Imat .- Imat')))
#this should be zero if the only issue is the order that integration is done.
#tis close for I, but not exact
display(Imat .- I2mat')
#%%
isapprox(Imat, Imat', rtol=1e-15)
isapprox(Wmat, Wmat', rtol=1e-15)

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


#%%

N = 10
A = zeros(ComplexF64, N, N)
B = zeros(ComplexF64, N, N)

for i in 2:N-1
    A[i, i] = 0.23271
    A[i, i+1] = 0.0581776 - 0.3im
    A[i+1, i] = 0.0581776 + 0.3im
    B[i, i] = 5.72958
    B[i, i+1] = -2.86479 - 0.3im
    B[i+1, i] = -2.86479 + 0.3im
end

A[1, 1] = 1.0
A[N, N] = 1.0
B[1, 1] = 1.0
B[N, N] = 1.0

A .- A'
B .- B'

ishermitian(A)
ishermitian(B)
isposdef(A)
isposdef(B)

eigen(A, B)

eigen(B)


A = [[2.0, 1+5im] [1-5im, 1.0]]
B = [[3.0, 2+5im] [2-5im, 1.0]]

ishermitian(A)
ishermitian(B)
isposdef(A)
#hmmm, so the actual eigenvalue problemo itself is not hermitian?
eigen(Hermitian(A), Hermitian(B))

eigen(Hermitian(A))
eigen(Hermitian(B))
isselfadjoint(A)

A .- adjoint(A)

isdef(A)
