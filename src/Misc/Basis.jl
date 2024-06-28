


"""
Arrays that distinguish the different Hermite basis functions, used for converting between the grid and the appropriate points indices in the matrices.
Defined here to match our construction, but used in Indexing.jl.
"""
const grid_id = [0, 0, 1, 1]
const basis_id = [0, 1, 0, 1]
const grid_id_r = [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1]
const grid_id_θ = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
const basis_id_r = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
const basis_id_θ = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]


"""
Creates the four Hermite basis functions, taken from wikipedia. Returns a 4xgp matrix for the 0th, 1st and 2nd derivative.

# Args
gp::Array{Float64} - Array of the Gaussian quadrature points.
"""
function hermite_basis(gp::Array{Float64})
    t = @. (gp + 1)/(2) #converts to correct range for spline
    H = zeros(4, length(gp))
    dH = zeros(4, length(gp))
    ddH = zeros(4, length(gp))


    ##############
    # node structure per element
    #element defined so locally x = [0, 1]
    # 0 -------- 1

    #each of the basis functions has either a value of 1 at one edge or a derivative of 1 at one edge
    #and then value of zero at the other three options.

    #ie first basis functions has H(0) = 1, H'(0) = 0, H(1) = 0, H'(1) = 0

    #node 1 H(0) = 1
    H[1, :] = @. 2*t^3 - 3*t^2 + 1

    #node 1 H'(0) = 1
    H[2, :] = @. 2*(t^3-2t^2 + t) #2 is from changin from -1, 1, to 1

    #node 2 H(1) = 0
    H[3, :] = @. -2t^3 + 3t^2

    #node 2 H'(1) = 0
    H[4, :] = @. 2*(t^3 - t^2) #2 is from changin from -1, 1, to 1

    #divide by 2 so each is 1 or 0 on the edges.
    dH[1, :] = @. (6*t^2 - 6*t) / 2
    dH[2, :] = @. 2*(3*t^2-4t + 1) / 2 #2 is from changin from -1, 1, to 1
    dH[3, :] = @. (-6t^2 + 6t) / 2
    dH[4, :] = @. 2*(3t^2 - 2t) / 2 #2 is from changin from -1, 1, to 1

    
    #these have to be divided again!
    ddH[1, :] = @. (12*t - 6) / 4
    ddH[2, :] = @. 2*(6*t-4) / 4 #2 is from changin from -1, 1, to 1
    ddH[3, :] = @. (-12t + 6) / 4
    ddH[4, :] = @. 2*(6t - 2) / 4 #2 is from changin from -1, 1, to 1

    return H, dH, ddH

end

"""
Creates the four Hermite basis functions, taken from wikipedia. Returns a 4xgp matrix for the 0th, 1st and 2nd derivative.

# Args
gp::Array{Float64} - Array of the Gaussian quadrature points.
"""
function hermite_basis(rgp::Array{Float64}, θgp::Array{Float64})

    Hr, dHr, ddHr = hermite_basis(rgp)
    Hθ, dHθ, ddHθ = hermite_basis(θgp)

    #probably could be more efficient, but this is only done once so shouldn't matter.
    S = combine_basis(Hr, Hθ, rgp, θgp)
    dSr = combine_basis(dHr, Hθ, rgp, θgp)
    dSθ = combine_basis(Hr, dHθ, rgp, θgp)
    ddSrr = combine_basis(ddHr, Hθ, rgp, θgp)
    ddSrθ = combine_basis(dHr, dHθ, rgp, θgp)
    ddSθθ = combine_basis(Hr, ddHθ, rgp, θgp)


    return S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ

end


#combines two 1d Hermite basis' into a 2d basis.
function combine_basis(Hr, Hθ, rgp, θgp)


    S = zeros(16, length(rgp), length(θgp))

    #recall from 1d
    ##############
    # node structure per element

    # 0 -------- 1

    #H[1] -> H(0) = 1
    #H[2] -> H'(0) = 1
    #H[3] -> H(1) = 1
    #H[4] -> H'(1) = 1


    ############### 
    #node labels for each element!

    # (0, 1) -------- (1, 1)
    #  |               |
    #  |               |
    #  |               |
    # (0, 0) -------- (1, 0)

    for i in 1:length(rgp), j in 1:length(θgp)

        # node (0, 0), S(0, 0) = 1
        S[1, i, j] = Hr[1, i] * Hθ[1, j]

        # Sr(0, 0) = 1 #stands for basis function that is derivative with respect to r.
        S[2, i, j] = Hr[2, i] * Hθ[1, j]

        #Sθ(0, 0) = 1
        S[3, i, j] = Hr[1, i] * Hθ[2, j]

        #Srθ(0, 0) = 1
        S[4, i, j] = Hr[2, i] * Hθ[2, j]

        #node (1, 0), S(1, 0) = 1
        S[5, i, j] = Hr[3, i] * Hθ[1, j]

        #Sr(1, 0) = 1
        S[6, i, j] = Hr[4, i] * Hθ[1, j]

        #Sθ(1, 0) = 1
        S[7, i, j] = Hr[3, i] * Hθ[2, j]

        #Srθ(1, 0) = 1
        S[8, i, j] = Hr[4, i] * Hθ[2, j]


        #node (0, 1), S(0, 1) = 1
        S[9, i, j] = Hr[1, i] * Hθ[3, j]

        #Sr(0, 1) = 1
        S[10, i, j] = Hr[2, i] * Hθ[3, j]

        #Sθ(0, 1) = 1
        S[11, i, j] = Hr[1, i] * Hθ[4, j]

        #Srθ(0, 1) = 1
        S[12, i, j] = Hr[2, i] * Hθ[4, j]


        #node (1, 1), S(1, 1) = 1
        S[13, i, j] = Hr[3, i] * Hθ[3, j]

        #Sr(1, 1) = 1
        S[14, i, j] = Hr[4, i] * Hθ[3, j]

        #Sθ(1, 1) = 1
        S[15, i, j] = Hr[3, i] * Hθ[4, j]

        #Srθ(1, 1) = 1
        S[16, i, j] = Hr[4, i] * Hθ[4, j]

    end

    return S


end




"""
Modifies the matrix storing the basis functions, Φ, to reflect the local derivatives being taken.

#Args
Φ::Array{ComplexF64, 3} - 9x4xgp matrix, stores the 9 derivtaives of the potential for each 4 Hermite basis functions at each gauss point.
H::Array{Float64, 2} - Hermite basis functions.
dH::Array{Float64, 2} - Derivtive of Hermite basis functions.
ddH::Array{Float64, 2} - Second derivtive of Hermite basis functions.
m::Int64 - Poloidal mode number, i.e. scale factor for derivative of Fourier basis.
n::Int64 - Toroidal mode number, i.e. scale factor for derivative of Fourier basis.
jac::Float64 - Jacobian of local to global transformation, not to be confused with geometric Jacobian.
"""
function create_local_basis!(Φ::Array{ComplexF64, 3}, H::Array{Float64, 2}, dH::Array{Float64, 2}, ddH::Array{Float64, 2}, m::Int64, n::Int64, jac::Float64)

    #Φ_s
    @. Φ[1, :, :] = dH / jac
    #Φ_θ
    @. Φ[2, :, :] = H * m * 1im
    #Φ_ζ
    @. Φ[3, :, :] = H * n * 1im
    #Φ_ss
    @. Φ[4, :, :] = ddH / jac^2
    #Φ_sθ
    @. Φ[5, :, :] = dH * m * 1im / jac
    #Φ_sζ   
    @. Φ[6, :, :] = dH * n * 1im / jac
    #Φ_θθ
    @. Φ[7, :, :] = H * (-m^2)
    #Φ_θζ
    @. Φ[8, :, :] = H * (-m*n)
    #Φ_ζζ
    @. Φ[9, :, :] = H * (-n^2)

    #test case with no derivative!
    #@. Φ[10, :, :] = H

end



#creates the locat basis for 2d fem case.
function create_local_basis!(Φ::Array{ComplexF64, 4}, S::Array{Float64, 3}, dSr::Array{Float64, 3}, dSθ::Array{Float64, 3}, ddSrr::Array{Float64, 3}, ddSrθ::Array{Float64, 3}, ddSθθ::Array{Float64, 3}, m::Int64, n::Int64, dr::Float64, dθ::Float64)

    rjac = dr / 2
    θjac = dθ / 2

    #unsure if we need to scale the non-deriv terms by e^(im) think most would cancel, but maybe not.
    #think they should cancel because every term will have e^imθ, regardless of derivative.
    #it will also be v tricky to scale by the theta coords.
    
    #Φ_s
    @. Φ[1, :, :, :] = dSr / rjac 
    #Φ_θ
    @. Φ[2, :, :, :] = dSθ / θjac + S * 1im * m
    #Φ_ζ
    @. Φ[3, :, :, :] = S * n * 1im 
    #Φ_ss
    @. Φ[4, :, :, :] = ddSrr / rjac^2 
    #Φ_sθ
    @. Φ[5, :, :, :] = ddSrθ / (rjac * θjac) + dSr * 1im * m / rjac
    #Φ_sζ   
    @. Φ[6, :, :, :] = dSr * n * 1im / rjac 
    #Φ_θθ
    @. Φ[7, :, :, :] = ddSθθ / θjac^2 + 2 * dSθ * 1im * m / θjac - m^2 * S
    #Φ_θζ
    @. Φ[8, :, :, :] = 1im * n * (dSθ / θjac + S * 1im * m)
    #Φ_ζζ
    @. Φ[9, :, :, :] = S * (-n^2) 

end

