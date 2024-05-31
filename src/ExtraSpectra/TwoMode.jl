"""
Function that computes the spectrum based on the reduced two mode TAE equation of Berk et al 1992.

# Args
rgrid::Array{Float64} - Radial grid
R::Float64 - Major radius
m::Int64 - First poloidal mode, considers modes m and m+1.
n::Int64 - Toroidal mode number.
δ::Float64 - Artificial resisitivity.
σ::Float64 - Target for shift-and-invert, typically TAE frequency.
"""
function two_mode(; rgrid::Array{Float64}, R0::Float64, m::Int64, n::Int64, δ::Float64, σ::Float64)

    m1 = m+1
    N = length(rgrid)

    gp = 5
    ξ, wg = gausslegendre(gp)
    H, dH, ddH = hermite_basis(ξ)

    W = zeros(ComplexF64, 4*N, 4*N)
    I = zeros(ComplexF64, 4*N, 4*N)

    Em = zeros(3, 4, gp)
    Em1 = zeros(3, 4, gp)

    Em[1, :, :] = H
    Em[2, :, :] = dH
    Em[3, :, :] = ddH

    #probably not needed.
    Em1[1, :, :] = H
    Em1[2, :, :] = dH
    Em1[3, :, :] = ddH

    if abs(m)==1
        bc1 = [2, 2*N-1]
    else
        bc1 = [1, 2*N-1]
    end
    if abs(m1)==1
        bc2 = [2*N+2, 4*N-1]
    else
        bc2 = [2*N+1, 4*N-1]
    end

    bcinds = vcat(bc1, bc2)
    rowsW = Array{Int64}(undef, 0)
    colsW = Array{Int64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)
    Idata = Array{ComplexF64}(undef, 0)
    rowsI = Array{Int64}(undef, 0)
    colsI = Array{Int64}(undef, 0)
    

    for i in 1:N-1

        r, dr = local_to_global(i, ξ, collect(rgrid))

        #Berk's case,
        #=
        q = @. 1 + r^2

        dens = ones(length(r))
        ddens = zeros(length(r))
        =#
        

        #Bowden singular

        #think this will be hard-coded to cover the Bowden case.
        #as passing in a density derivative does not match other cases.
        
        q = @. 1.0 + 2.0*r^2

        dens = @. 1/2*(1-tanh((r-0.7)/0.05))

        ddens = @. -10 * sech((r-0.07)/0.05)^2
        

        #Axel 2012
        #=
        q = @. 1.05 + 0.55*r^2

        #dens = ones(length(r)) #for testing.
        #ddens = zeros(length(r))
        dens = @. 1/2*(1-tanh((r-0.8)/0.1))
        ddens = @. -5 * sech(8-10*r)^2
        =#


        km = @. 1/R0*(m/q + n)

        #kmp = @. -m*qp/(R*q^2)

        #kmpp = @. m * (2*qp^2 - q*qpp)/(R*q^3)

        km1 = @. 1/R0*(m1/q + n)

        #km1p = @. -m1*qp/(R*q^2)

        #km1pp = @. m1 * (2*qp^2 - q*qpp)/(R*q^3)

        jac = dr/2

        ϵ = @. 5*r/(2*R0)
        #Fu Dam case 
        #this gives the exact frequency that Berk quotes in fig 1...
        #ϵ = @. 3*r/(2*R)

        #cylindrical limit!
        #ϵ = zeros(length(r))

        


        #equation 1
        for trial in 1:4

            trial_ind_m = grid_to_ind_coupled(i, N, trial, 1)
            trial_ind_m1 = grid_to_ind_coupled(i, N, trial, 2)

            for test in 1:4
                test_ind_m = grid_to_ind_coupled(i, N, test, 1)
                test_ind_m1 = grid_to_ind_coupled(i, N, test , 2)

                for j in 1:gp
                    #Eq 1 I
                    push!(rowsI, test_ind_m)
                    push!(colsI, trial_ind_m)
                    push!(Idata, jac * Em[2, test, j] * r[j]^3 * Em[2, trial, j] * 1/jac^2 * wg[j] * dens[j])

                    push!(rowsI, test_ind_m)
                    push!(colsI, trial_ind_m)
                    push!(Idata, jac * (m^2-1) * Em[1, test, j] * r[j] * Em[1, trial, j] * wg[j] * dens[j])

                    push!(rowsI, test_ind_m)
                    push!(colsI, trial_ind_m)
                    push!(Idata, -jac * Em[1, test, j] * r[j]^2 * ddens[j] * Em[1, trial, j] * wg[j])

                    push!(rowsI, test_ind_m)
                    push!(colsI, trial_ind_m1)
                    push!(Idata, jac * Em[2, test, j] * r[j]^3 * ϵ[j] * Em1[2, trial, j] * 1/jac^2 * wg[j] * dens[j])

                    push!(rowsI, test_ind_m)
                    push!(colsI, trial_ind_m)
                    push!(Idata, -jac * 1im * δ * (Em[3, test, j] / jac^2 + 1/r[j] * Em[2, test, j] / jac - m^2 / r[j]^2 * Em[1, test, j]) * dens[j] * (r[j] * Em[3, trial, j] / jac^2 + 3 * Em[2, trial, j] / jac + Em[1, trial, j] /r[j] - m^2/r[j] * Em[1, trial, j]))

                    #eq1 W
                    push!(rowsW, test_ind_m)
                    push!(colsW, trial_ind_m)
                    push!(Wdata, jac * Em[2, test, j] * r[j]^3 * km[j]^2 * Em[2, trial, j] * 1/jac^2 * wg[j])

                    push!(rowsW, test_ind_m)
                    push!(colsW, trial_ind_m)
                    push!(Wdata, jac * (m^2-1) * Em[1, test, j] * r[j] * km[j]^2 * Em[1, trial, j] * wg[j]) 
                    

                    #eq2 I
                    push!(rowsI, test_ind_m1)
                    push!(colsI, trial_ind_m1)
                    push!(Idata, jac * Em1[2, test, j] * r[j]^3 * Em1[2, trial, j] * 1/jac^2 * wg[j] * dens[j])

                    push!(rowsI, test_ind_m1)
                    push!(colsI, trial_ind_m1)
                    push!(Idata, jac * (m1^2-1) * Em1[1, test, j] * r[j] * Em1[1, trial, j] * wg[j] * dens[j])

                    push!(rowsI, test_ind_m1)
                    push!(colsI, trial_ind_m1)
                    push!(Idata, -jac * Em1[1, test, j] * r[j]^2 * ddens[j] * Em1[1, trial, j] * wg[j])

                    push!(rowsI, test_ind_m1)
                    push!(colsI, trial_ind_m)
                    push!(Idata, jac * Em1[2, test, j] * r[j]^3 * ϵ[j] * Em[2, trial, j] * 1/jac^2 * wg[j] * dens[j])

                    push!(rowsI, test_ind_m1)
                    push!(colsI, trial_ind_m1)
                    push!(Idata, -jac * 1im * δ * (Em1[3, test, j] / jac^2 + 1/r[j] * Em1[2, test, j] / jac - m1^2 / r[j]^2 * Em1[1, test, j]) * dens[j] * (r[j] * Em1[3, trial, j] / jac^2 + 3 * Em1[2, trial, j] / jac + Em1[1, trial, j] /r[j] - m1^2/r[j] * Em1[1, trial, j]))



                    #eq2 W
                    push!(rowsW, test_ind_m1)
                    push!(colsW, trial_ind_m1)
                    push!(Wdata, jac * Em1[2, test, j] * r[j]^3 * km1[j]^2 * Em1[2, trial, j] * 1/jac^2 * wg[j])

                    push!(rowsW, test_ind_m1)
                    push!(colsW, trial_ind_m1)
                    push!(Wdata, jac * (m1^2-1) * Em1[1, test, j] * r[j] * km1[j]^2 * Em1[1, trial, j] * wg[j]) 


                end
            end
        end
    end

    #this may be slow af, but should work as an awful solution...
    #hack solution to keep the matrices sparse, works, but probably not very well.
    rowsbcI = findall(x->x in bcinds, rowsI)
    colsbcI = findall(x->x in bcinds, colsI)
    Idata[rowsbcI] .= 0.0
    Idata[colsbcI] .= 0.0
    rowsbcW = findall(x->x in bcinds, rowsW)
    colsbcW = findall(x->x in bcinds, colsW)
    Wdata[rowsbcW] .= 0.0
    Wdata[colsbcW] .= 0.0
    for i in bcinds
        push!(rowsI, i)
        push!(colsI, i)
        push!(rowsW, i)
        push!(colsW, i)
        push!(Wdata, 1.0)
        push!(Idata, 1.0)
    end
    #this is fast as fk boi, now we have sparse fellas.
    W = sparse(rowsW, colsW, Wdata)
    I = sparse(rowsI, colsI, Idata)
    
    

    #display(Matrix(I))
    #display(I)
    #even with current term, matrix seems `pretty` Hermitian, ie ~e-15 but now we get negative evals???
    #display(maximum(W - W')) #so these matrices are Hermitian to a tolerance of ~e-15
    #ω2, evals = eigen(Matrix(Hermitian(W)), Matrix(Hermitian(I)))
    ω2, evals = eigs(W, I, nev=5, ritzvec=true, sigma=σ)
    #ω2, evals = eigs(W, I, nev=5, ritzvec=true, sigma=(0.39/R)^2)

    Emsol = evals[1:2:2*N, :]
    Em1sol = evals[2*N+1:2:end, :]
    return R0 .*sqrt.(ω2), Emsol, Em1sol

end


"""
Maps the grid indexes to matrix indexes for the two mode case.

# Args
i::Int64 - Grid index
N::Int64 - Number of grid points
h::Int64 - Basis function index
eqn::Int64 - Int labelling which of the two equations.
"""
function grid_to_ind_coupled(i::Int64, N::Int64, h::Int64, eqn::Int64)

    grid_id = [0, 0, 1, 1]

    basis_id = [0, 1, 0, 1]

    return 2 * (eqn-1)*N + 1 + 2*(i-1) + basis_id[h] + 2*grid_id[h]

end

"""
Crudely reconstructs the continuum based on eigenfunctions for two mode case.

# Args
ω::Array{ComplexF64} - Eigenvalues
rgrid::Array{Float64} - Radial grid
Em::Array{ComplexF64, 2} - First eigenfunction
Em1::Array{ComplexF64, 2} - Second eigenfunction
"""
function Em_to_cont(; ω::Array{ComplexF64}, rgrid::Array{Float64}, Em::Array{ComplexF64, 2}, Em1::Array{ComplexF64, 2})
    rdata = zeros(length(ω));
    omdata = zeros(ComplexF64, length(ω));
    col = zeros(length(ω));
    for i in 1:1:length(ω)

        rm = argmax(abs.(rgrid .*Em[:, i]))
        rm1 = argmax(abs.(rgrid .*Em1[:, i]))
        if abs(Em[rm, i]) > abs(Em1[rm1, i])

            rdata[i] = rgrid[rm]
            col[i] = 0.0
        else
            rdata[i] = rgrid[rm1]
            col[i] = 1.0
        end

        omdata[i] = ω[i]
    end
    return rdata, omdata, col
end