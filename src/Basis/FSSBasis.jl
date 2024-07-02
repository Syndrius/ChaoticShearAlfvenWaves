
"""
    hermite_basis(gp::Array{Float64})

Creates the four Hermite basis functions, taken from wikipedia. Returns a 4xgp matrix for the 0th, 1st and 2nd derivative.
"""
function hermite_basis(gp::Array{Float64})
    t = @. (gp + 1)/(2) #converts to correct range
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
    create_local_basis!(Φ::Array{ComplexF64, 3}, H::Array{Float64, 2}, dH::Array{Float64, 2}, ddH::Array{Float64, 2}, m::Int64, n::Int64, jac::Float64)

Modifies the matrix storing the basis functions, Φ, to reflect the local derivatives being taken.

### Args
Φ::Array{ComplexF64, 3} - 9x4xgp matrix, stores the 9 derivtaives of the potential for each 4 Hermite basis functions at each gauss point.
H::Array{Float64, 2} - Hermite basis functions.
dH::Array{Float64, 2} - Derivtive of Hermite basis functions.
ddH::Array{Float64, 2} - Second derivtive of Hermite basis functions.
m::Int64 - Poloidal mode number, i.e. scale factor for derivative of Fourier basis.
n::Int64 - Toroidal mode number, i.e. scale factor for derivative of Fourier basis.
jac::Float64 - Jacobian of local to global transformation, not to be confused with geometric Jacobian.
"""
function create_local_basis!(Φ::Array{ComplexF64, 3}, H::Array{Float64, 2}, dH::Array{Float64, 2}, ddH::Array{Float64, 2}, m::Int64, n::Int64, jac::Float64)

    #the middle index is chosen to increase efficiency with numerical integration.

    #Φ_s
    @. Φ[:, 1, :] = dH / jac
    #Φ_θ
    @. Φ[:, 2, :] = H * m * 1im
    #Φ_ζ
    @. Φ[:, 3, :] = H * n * 1im
    #Φ_ss
    @. Φ[:, 4, :] = ddH / jac^2
    #Φ_sθ
    @. Φ[:, 5, :] = dH * m * 1im / jac
    #Φ_sζ   
    @. Φ[:, 6, :] = dH * n * 1im / jac
    #Φ_θθ
    @. Φ[:, 7, :] = H * (-m^2)
    #Φ_θζ
    @. Φ[:, 8, :] = H * (-m*n)
    #Φ_ζζ
    @. Φ[:, 9, :] = H * (-n^2)

end
