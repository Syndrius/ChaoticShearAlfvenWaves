using Bessels
using FunctionZeros
#%%
function anal_sol(i, m, n, xgrid, ygrid, zgrid)
    #ith zero
    #mth bessel function and poloidal mode number
    #n toroidal mode number

    Nx = length(xgrid)
    Ny = length(ygrid)
    Nz = length(zgrid)

    sol = zeros(Nx, Ny, Nz)

    j0 = besselj_zero(m, i)

    for (j, x) in enumerate(xgrid), (k, y) in enumerate(ygrid), (l, z) in enumerate(zgrid)

        #unsure about the sin and cos.
        sol[j, k, l] = besselj(m, j0*x)*(sin(m*y) + cos(m*y))*(sin(n*z) + cos(n*z))
    end

    return sol
end

function anal_evals(irange, mrange, nrange)

    evals = Float64[]
    code = []
    for i in irange, m in mrange
        j0 = besselj_zero(m, i)

        for n in nrange
            push!(evals, j0^2 + n^2)
            push!(code, (i, m, n))
        end
    end

    perm = sortperm(evals)

    return evals[perm], code[perm]
end

function anal_eval(i, m, n)
    j0 = besselj_zero(m, i)

    return j0^2 + n^2

end

