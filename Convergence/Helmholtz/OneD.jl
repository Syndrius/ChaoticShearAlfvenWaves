
#test case solving the Helmholtz equation

function compute_linear_bcs(N::Int64)
    #return [N]
    return [1, N]
end

function compute_hermite_bcs(N::Int64)
    return [1, 2*N-1]
end

function cubic_mat_size(N::Int64)
    return 2*N
end

#assume 1d for now.
function linear(N::Int64, gp::Int64)

    xgrid = LinRange(0, π, N)
    xgrid = collect(LinRange(0, 1, N))

    ξ, wg = gausslegendre(gp)

    S = linear_basis(ξ)

    M = zeros(Float64, N, N)

    A = zeros(Float64, N, N)

    bcs = compute_linear_bcs(N)
    n = 1

    for i in 1:N-1

        x, dx = linear_local_to_global(i, ξ, xgrid)

        jac = dx / 2

        for trialx in 1:2

            r_ind = grid_to_index_linear(i, trialx, N)

            if r_ind in bcs
                continue
            end

            for testx in 1:2

                l_ind = grid_to_index_linear(i, testx, N)

                if l_ind in bcs
                    continue
                end

                for kx in 1:gp
                    #Cartesian Helmholtz
                    #M[l_ind, r_ind] += jac * (S.dH[testx, kx] * S.dH[trialx, kx]) / jac^2 * wg[kx] 
                    #A[l_ind, r_ind] += jac * S.H[testx, kx] * S.H[trialx, kx] * wg[kx]

                    #radial part of 2d equaiton.
                    M[l_ind, r_ind] += jac * (-S.dH[testx, kx] * S.dH[trialx, kx] / jac^2
                                              + S.H[testx, kx] / x[kx] * S.dH[trialx, kx] / jac
                                              - n^2 * S.H[testx, kx] * S.H[trialx, kx] / x[kx]^2) * wg[kx]

                    A[l_ind, r_ind] += -jac * (S.H[testx, kx] * S.H[trialx, kx]) * wg[kx]

                end
            end
        end
    end

    #for i in 1:N
    #    A[i, i] /= 2
    #end

    for i in bcs
        M[i, i] = 1.0
        A[i, i] = 1.0
    end

    #display(M)
    #display(A)

    return eigen(M, A)
    return eigen(Hermitian(M), Hermitian(A))
end#
#
function h00(t::Float64)

    return 2*t^3 - 3*t^2 + 1
end

function h10(t::Float64)
    #return 2 * (t^3-2*t^2+t)
    return t*(1-t)^2
end

function h01(t::Float64)
    return -2t^3+3t^2
    #return t^2*(3-2*t)
end

function h11(t::Float64)
    #return 2*(t^3-t^2)
    return t^2*(t-1)
end

struct NHB1d
    H :: Array{Float64, 2}
    dH :: Array{Float64, 2}
    ddH :: Array{Float64, 2}
end

function dh00(t)

    return 6*t^2 - 6*t
end

function dh10(t)
    return (1-t)^2 - 2*t*(1-t)
end

function dh01(t)
    return -6*t^2 + 6*t
end

function dh11(t)
    return 2*t*(t-1) + t^2
end

function ddh00(t)
    return 12*t - 6
end

function ddh10(t)
    return -4*(1-t)+2*t
end

function ddh01(t)
    return -12*t + 6
end

function ddh11(t)
    return 4*t + 2*(t-1)
end



function new_hermite_basis(ξ::Array{Float64})

    S = NHB1d(zeros(4, length(ξ)), zeros(4, length(ξ)), zeros(4, length(ξ)))

    #our code looks to define this for ξ, so they are defined from -1, to 1.
    #unsure if that is better! -> perhaps it is, then we only need a single local grid
    #i.e. from -1 to 1.
    #also unsure why we are dividing by 2.
    #seems like that is not needed and doesn't really make any sense.
    #that would explain the dx/2 needed as the extra factor.

    #this is essentially identical now, but this should be clearer.
    #the main feature is that we have shifted the basis function to be defined across -1, 1
    #to match the domain for gaussian quadrature, hence our local domain is always -1 to 1.
    t = @. (ξ + 1) / 2
    dtdξ = 1 / 2
    dξdt = 2
    S.H[1, :] = h00.(t)
    S.H[2, :] = h10.(t) * dξdt #reflects change of coords to ξ
    S.H[3, :] = h01.(t)
    S.H[4, :] = h11.(t) * dξdt #reflects change of coords to ξ

    S.dH[1, :] = dh00.(t) * dtdξ
    S.dH[2, :] = dh10.(t) #* dξdt * dtdξ #no net result.
    S.dH[3, :] = dh01.(t) * dtdξ
    S.dH[4, :] = dh11.(t) #* dξdt * dtdξ #no net result.

    S.ddH[1, :] = ddh00.(t) * dtdξ * dtdξ
    S.ddH[2, :] = ddh10.(t) * dtdξ #* dξdt * dtdξ #no net result.
    S.ddH[3, :] = ddh01.(t) * dtdξ * dtdξ
    S.ddH[4, :] = ddh11.(t) * dtdξ #* dξdt * dtdξ #no net result.

   return S

end

function cubic(N::Int64, gp::Int64, n=1)


    xgrid = collect(LinRange(0, π, N))
    #xgrid = collect(LinRange(0, 1, N))

    ξ, wg = gausslegendre(gp)

    S = new_hermite_basis(ξ)
    S2 = hermite_basis(ξ)

    display(S.ddH .- S2.ddH)

    mat_size = cubic_mat_size(N)

    M = zeros(Float64, mat_size, mat_size)

    A = zeros(Float64, mat_size, mat_size)

    bcs = compute_hermite_bcs(N)

    for i in 1:N-1

        #this is actually the same as full case.
        x, dx = linear_local_to_global(i, ξ, xgrid)

        jac = dx / 2

        #this is an awful solution but may be easiest.
        #ideally we an put this into local_basis I think! -> should be plausible, just need to create 
        #a lil matrix to add into the local basis
        #we need to understand why it takes this form though!
        #dx/2 seems to work, but difficult to tell if this is ok
        #need to understand why this is.
        #dx/2 is because we are shifting from [-1, 1] domain to [x_i, x_{i+1}] domain
        #so the tangent shape functions need an extra scaling of dx/2 from the chain rule.
        #this is not the way to put this into the function!
        #should be done in the local basis step./
        sf = [1.0, dx/2, 1.0, dx/2]

        for trialx in 1:4

            r_ind = grid_to_index_cubic(i, trialx, N)

            if r_ind in bcs
                continue
            end

            for testx in 1:4

                l_ind = grid_to_index_cubic(i, testx, N)

                if l_ind in bcs
                    continue
                end

                for kx in 1:gp

                    #looks like we will be able to artifically shift to higher derivs with no problem.
                    #M[l_ind, r_ind] += jac * (S.ddH[testx, kx] * S.ddH[trialx, kx] / jac^4) * wg[kx]
                    #A[l_ind, r_ind] += jac * (S.dH[testx, kx] * S.dH[trialx, kx]) / jac^2 * wg[kx] 

                    #normal cartesian case
                    M[l_ind, r_ind] += jac * (S.dH[testx, kx] * S.dH[trialx, kx]) / jac^2 * wg[kx] * sf[testx] * sf[trialx]
                    A[l_ind, r_ind] += jac * S.H[testx, kx] * S.H[trialx, kx] * wg[kx] * sf[testx] * sf[trialx]

                    #polar cas with a given n.
                    #M[l_ind, r_ind] += jac * (-S.dH[testx, kx] * S.dH[trialx, kx] / jac^2
                    #                          + S.H[testx, kx] / x[kx] * S.dH[trialx, kx] / jac
                    #                          - n^2 * S.H[testx, kx] * S.H[trialx, kx] / x[kx]^2) * wg[kx]

                    #A[l_ind, r_ind] += -jac * (S.H[testx, kx] * S.H[trialx, kx]) * wg[kx]

                end
            end
        end
    end

    #ok for some reason this makes the result much better
    #looks like we do have a serious problemo!
    #unsure if the result is actually corret here though!
    #this does cook the eigenvalues tho...
    #oh dear, something big is wrong...


    display(A)
    for i in bcs
        M[i, i] = 1.0
        A[i, i] = 1.0
    end


    return eigen(M, A)
    return eigen(Hermitian(M), Hermitian(A))

end

#so this is way better than the fem approach, eigenvalues are a little less acurate, 
#but similar to linear, but the eigenfunctions are bang on, including the amplitude
#which is not true for fem.
function basic(N::Int64)
    #just fill matrices in directly to see if the amplitude is ok
    A = zeros(N, N)
    B = zeros(N, N)

    #dx = 1.0 / (N-1)
    dx = π / (N-1)

    for i in 1:N

        if (i==1)
            A[1, 1] = -2.0 / dx^2
            A[1, 2] = 1.0 / dx^2
        elseif (i==N)
            A[N, N-1] = 1.0 / dx^2
            A[N, N] = - 2.0 / dx^2
        else
            A[i, i-1] = 1.0 / dx^2
            A[i, i] = -2.0 / dx^2
            A[i, i+1] = 1.0 / dx^2
        end
        B[i, i] = -1.0
    end

    #boundaries.
    A[1, 1] = 1.0
    A[1, 2] = 0.0

    B[1, 1] = 0.0

    A[N, N] = 1.0
    A[N, N-1] = 0.0

    B[N, N] = 0.0

    #display(A)
    return eigen(A, B)
end

function direct(N::Int64)
    A = zeros(N, N)
    B = zeros(N, N)

    xgrid = LinRange(0, π, N)

    l = xgrid[2] - xgrid[1] #const for uniform grid

    for i in 2:N-1

        A[i, i] = 2.0 / l
        B[i, i] = 4.0 * l / 6
        A[i, i+1] = -1.0 / l
        B[i, i+1] = 1.0 * l / 6
        A[i, i-1] = -1.0 / l
        B[i, i-1] = 1.0 * l / 6
    end

    A[1, :] .= 0.0
    A[:, 1] .= 0.0
    A[1, 1] = 1.0
    A[N, :] .= 0.0
    A[:, N] .= 0.0
    A[N, N] = 1.0
    B[1, :] .= 0.0
    B[:, 1] .= 0.0
    B[1, 1] = 1.0
    B[N, :] .= 0.0
    B[:, N] .= 0.0
    B[N, N] = 1.0

    #display(A)
    #display(B)

    return eigen(A, B)

end


#this may be the best place to understand our aproach
#as this skips the local to global transformation
#meaning we can focus purely on the basis functions!
function no_integration(N::Int64, gp::Int64, n=1)

    xgrid = collect(LinRange(0, π, N))
    #xgrid = collect(LinRange(0, 1, N))

    dx = xgrid[2] - xgrid[1] #const for this case.

    #ξ, wg = gausslegendre(gp)

    #S = new_hermite_basis(ξ)
    #S = hermite_basis(ξ)

    mat_size = cubic_mat_size(N)

    M = zeros(Float64, mat_size, mat_size)

    A = zeros(Float64, mat_size, mat_size)

    bcs = compute_hermite_bcs(N)

    Ist = [[13/35, 9/70, 11/210, -13/420]  [9/70, 13/35, 13/420, -11/210]  [11/210, 13/420, 1/105, -1/140]  [-13/420, -11/210, -1/140, 1/105]]
    #Ist = [[13/35+ 9/70+ 11/210-13/420, 9/70+ 13/35+ 13/420+ -11/210]  [11/210+ 13/420+ 1/105+ -1/140, -13/420+ -11/210+ -1/140+ 1/105]]
    Wst = [[6/5, -6/5, 1/10, 1/10]  [-6/5, 6/5, -1/10, -1/10]  [1/10, -1/10, 2/15, -1/30]  [1/10, -1/10, -1/30, 2/15]]
    Ist = [[13/35, 11/210, 9/70, -13/420]  [9/70, 13/420, 13/35, -11/210]  [11/210, 1/105, 13/420, -1/140]  [-13/420, -1/140, -11/210, 1/105]]
    Wst = [[6/5, 1/10, -6/5, 1/10]  [-6/5, -1/10, 6/5, -1/10]  [1/10, 2/15, -1/10, -1/30]  [1/10, -1/30, -1/10, 2/15]]
    Wst = [[3/10, -3/10, 1/40, 1/40]  [-3/10, 3/10, -1/40, -1/40]  [1/40, -1/40, 1/30, -1/120]  [1/40, -1/40, -1/120, 1/30]]

    Ist = [[13/35, 11/210, 9/70, -13/420]  [11/210, 1/105, 13/420, -1/140]  [9/70, 13/420, 13/35, -11/210]  [-13/420, -1/140, -11/210, 1/105]]
    Wst = [[3/10, 1/40, -3/10, 1/40]  [1/40, 1/30, -1/40, -1/120]  [-3/10, -1/40, 3/10, -1/40]  [1/40, -1/120, -1/40, 1/30]]

    #this seems to be a very important step for getting the correct ratio
    #I imagine this will be a bit more important where dx is not constant.
    #we probbaly need to test that!
    Ist = [[13/35, 11/210*dx, 9/70, -13/420*dx]  [11/210*dx, 1/105*dx^2, 13/420*dx, -1/140*dx^2]  [9/70, 13/420*dx, 13/35, -11/210*dx]  [-13/420*dx, -1/140*dx^2, -11/210*dx, 1/105*dx^2]]
    Wst = [[3/10, 1/40*dx, -3/10, 1/40*dx]  [1/40*dx, 1/30*dx^2, -1/40*dx, -1/120*dx^2]  [-3/10, -1/40*dx, 3/10, -1/40*dx]  [1/40*dx, -1/120*dx^2, -1/40*dx, 1/30*dx^2]]
    #=
    Wst[1, 2] /= dx / 2
    Wst[2, 1] /= dx / 2
    Wst[2, 2] /= dx^2 / 4
    Wst[3, 2] /= dx / 2
    Wst[2, 3] /= dx / 2
    Wst[2, 4] /= dx^2 / 4
    Wst[4, 2] /= dx^2 / 4
    Wst[4, 4] /= dx^2 / 4
    Wst[1, 4] /= dx / 2
    Wst[4, 1] /= dx / 2
    Wst[3, 4] /= dx / 2
    Wst[4, 3] /= dx / 2

    Ist[1, 2] /= dx / 2
    Ist[2, 1] /= dx / 2
    Ist[2, 2] /= dx^2 / 4
    Ist[3, 2] /= dx / 2
    Ist[2, 3] /= dx / 2
    Ist[2, 4] /= dx^2 / 4
    Ist[4, 2] /= dx^2 / 4
    Ist[4, 4] /= dx^2 / 4
    Ist[1, 4] /= dx / 2
    Ist[4, 1] /= dx / 2
    Ist[3, 4] /= dx / 2
    Ist[4, 3] /= dx / 2
    =#

    for i in 1:N-1

        #this is actually the same as full case.
        #x, dx = linear_local_to_global(i, ξ, xgrid)

        #jac = dx / 2

        li = 2*(i-1)+1
        ri = li+3
        display(li)
        display(ri)
        

        #display(A[li:ri, li:ri])

        #M[2*(i-1)+1:2*i+1, 2*(i-1)+1:2*i+1] .+= Wst
        A[li:ri, li:ri] .+= Ist 
        M[li:ri, li:ri] .+= Wst .* 4 / dx^2


    end
    display(A)
    #M[1:4, 1:4] .+= Ist
    A[1, :] .= 0.0
    A[:, 1] .= 0.0
    A[1, 1] = 1.0
    M[1, :] .= 0.0
    M[:, 1] .= 0.0
    M[1, 1] = 1.0

    A[2*N-1, :] .= 0.0
    A[:, 2*N-1] .= 0.0
    A[2*N-1, 2*N-1] = 1.0
    M[2*N-1, :] .= 0.0
    M[:, 2*N-1] .= 0.0
    M[2*N-1, 2*N-1] = 1.0

    return eigen(M, A)
    
end
