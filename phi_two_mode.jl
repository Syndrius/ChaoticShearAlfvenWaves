

#direct implementation of Axel's equation for phi
#hopefully we can replicated Axel's results
#may not worry about doing this, seems a bit hard tbh!

using FastGaussQuadrature
using MID
using LinearAlgebra
using Arpack


function grid_to_ind_coupled(i, N, h, eqn)

    grid_id = [0, 0, 1, 1]

    basis_id = [0, 1, 0, 1]

    return 2 * (eqn-1)*N + 1 + 2*(i-1) + basis_id[h] + 2*grid_id[h]

end

#TODO!
function Em_to_cont(ω, rgrid, Em, Em1)
    rdata = zeros(length(ω));
    omdata = zeros(ComplexF64, length(ω));
    col = zeros(length(ω));
    for i in 1:1:length(ω)
        #this is much much better for cylinder case.
        #maybe this would be better in general then??
        rm = argmax(abs.(rgrid .*Em[:, i]))
        rm1 = argmax(abs.(rgrid .*Em1[:, i]))
        #rm = argmax(abs.(Em[:, i]))
        #rm1 = argmax(abs.(Em1[:, i]))
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


function phi_2mode(N, m)
    #I think unfortunatly we have to do this now, think we are fked until we can replicated Axel's code.
    rgrid = clustered_grid(N, 0.91, 0.98, 0.25) #will just assume Axel's case for this.
    gp = 5
    ξ, wg = gausslegendre(gp)
    H, dH, ddH = MID.Misc.hermite_basis(ξ)

    W = zeros(ComplexF64, 4*N, 4*N)
    I = zeros(ComplexF64, 4*N, 4*N)

    ϕm = zeros(3, 4, gp)
    ϕm1 = zeros(3, 4, gp)

    ϕm[1, :, :] = H
    ϕm[2, :, :] = dH
    ϕm[3, :, :] = ddH

    #probably not needed.
    ϕm1[1, :, :] = H
    ϕm1[2, :, :] = dH
    ϕm1[3, :, :] = ddH

    if abs(m)==0
        bc1 = [2, 2*N-1]
    else
        bc1 = [1, 2*N-1]
    end
    if abs(m1)==0
        bc2 = [2*N+2, 4*N-1]
    else
        bc2 = [2*N+1, 4*N-1]
    end

    bcinds = vcat(bc1, bc2)

    rows = Array{Int64}(undef, 0)
    cols = Array{Int64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)
    Idata = Array{ComplexF64}(undef, 0)

    for i in 1:N-1

        r, dr = MID.Misc.local_to_global(i, ξ, collect(rgrid))

        Δp = r/(4*R0)

        ϵg = 2*Δp

        ϵb = - r/R0 

        ϵ = ϵg/2 - 2 * ϵb # 5/2 * r/R good sign.

        q = @. 1.05 + 0.55*r^2

        dens = @. 1/2*(1-tanh((r-0.8)/0.1))

        #think we can do the weak form so that we don't need the derivative of km...
        km = (m/q + n)/R0 #will need to confirm this? may be -n instead
        km1 = (m1/q + n)/R0

        for trial in 1:4

            trial_ind_m = grid_to_ind_coupled(i, N, trial, 1)
            trial_ind_m1 = grid_to_ind_coupled(i, N, trial, 2)
            for test in 1:4
                test_ind_m = grid_to_ind_coupled(i, N, test, 1)
                test_ind_m1 = grid_to_ind_coupled(i, N, test, 2)

                #I will be left hand side of Axel's equation.
                #weak form will defo multiply by r.
                #we will probably have to derive Axels result from equation given in paper.
                Isum11 = 0

                Isum11 += #TODO

                Wsum11 = 0

                Wsum11 += #TODO


                push!(rows, test_ind_m)
                push!(cols, trial_ind_m)
                push!(Idata, Isum11)
                push!(Wdata, Wsum11)


                Isum12 = 0

                Isum12 += #TODO

                Wsum12 = 0

                Wsum12 += #TODO

                push!(rows, test_ind_m)
                push!(cols, trial_ind_m1)
                push!(Idata, Isum12)
                push!(Wdata, Wsum12)


                Isum21 = 0

                Isum21 += #TODO

                Wsum21 = 0

                Wsum21 += #TODO

                push!(rows, test_ind_m1)
                push!(cols, trial_ind_m)
                push!(Idata, Isum21)
                push!(Wdata, Wsum21)



                Isum22 = 0

                Isum22 += #TODO

                Wsum22 = 0

                Wsum22 += #TODO

                push!(rows, test_ind_m1)
                push!(cols, trial_ind_m1)
                push!(Idata, Isum22)
                push!(Wdata, Wsum22)

            end
        end
    end

    #still an awful way to do bc's but seems to work
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


    #then solve...

        
end