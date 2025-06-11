#are the matrices actually correct.
#seems like to within machine precision, they are fine, but maybe that is where the problemos are occuring?
#perhaps the small numerical noise is causing issues?
#not sure that is a fixable problemo
using MID
using LinearAlgebra

#%%
#first define the problem
geo = init_geo(R0=4.0)

#rgrid = collect(LinRange(0, 1, N));
prob = init_problem(q=fu_dam_q, geo=geo); 
#%%
#then create the grids
Nr = 100;
rgrid = init_grid(type=:rf, N=Nr)
θgrid = init_grid(type=:as, N = 2, start = 1)
ζgrid = init_grid(type=:as, N = 1, start = -1)
grids = init_grids(rgrid, θgrid, ζgrid);

W, I = MID.Construct.construct(prob, grids)

#%%

ishermitian(W)
isadjoint(W)
isposdef(W)
issymmetric(W)

isapprox(W, W', rtol=1e-15)

display(Matrix(W)[1:5, 1:5])
display(Matrix(W')[1:5, 1:5])

eigvals(Matrix(W))
