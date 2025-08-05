using Bessels
using FunctionZeros
#%%

function rad_anal_sol(i, m, xgrid)
    #ith zero
    #mth bessel function and poloidal mode number
    #n toroidal mode number

    Nx = length(xgrid)

    sol = zeros(Nx)

    j0 = besselj_zero(m, i)

    for (j, x) in enumerate(xgrid)

        #unsure about the sin and cos.
        sol[j] = besselj(m, j0*x)
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

function zero_crossings(f, N) 
    nflip = 0
    for i in 2:N-2 #ignore start and finish as we know they are zeros
        if sign(real(f[i])) ≠ sign(real(f[i+1])) #obvs will need the m, n combo we are using! -> annoying for sign of m, n
            nflip += 1
        end
    end
    #plus 1 here is to match bessel definition.
    #i..e the end point is considered a zero.
    return nflip + 1
end
