
"""
Struct for the hermite basis functions, for 1d finite elements.

### Fields
- H::Array{Float64, 2} - Shape functions, size [4, gp]
- dH::Array{Float64, 2} - Derivatives, size [4, gp]
- ddH::Array{Float64, 2} - Second derivatives, size [4, gp]
"""
struct HB1d
    H :: Array{Float64, 2}
    dH :: Array{Float64, 2}
    ddH :: Array{Float64, 2}
end


"""
    hermite_basis(gp::Array{Float64})

Creates the four Hermite basis functions, taken from wikipedia. Returns a 4xgp matrix for the 0th, 1st and 2nd derivative.
"""
function hermite_basis(gp::Array{Float64})
    t = @. (gp + 1)/(2) #converts to correct range

    #creates a struct to store the hermiete shape functions, their derivative and second derivative.
    S = HB1d(zeros(4, length(gp)), zeros(4, length(gp)), zeros(4, length(gp)))

    ##############
    # node structure per element
    #element defined so locally x = [0, 1]
    # 0 -------- 1

    #each of the basis functions has either a value of 1 at one edge or a derivative of 1 at one edge
    #and then value of zero at the other three options.

    #ie first basis functions has H(0) = 1, H'(0) = 0, H(1) = 0, H'(1) = 0

    #node 1 H(0) = 1
    S.H[1, :] = @. 2*t^3 - 3*t^2 + 1

    #node 1 H'(0) = 1
    S.H[2, :] = @. 2*(t^3-2t^2 + t) #2 is from changin from -1, 1, to 1

    #node 2 H(1) = 0
    S.H[3, :] = @. -2t^3 + 3t^2

    #node 2 H'(1) = 0
    S.H[4, :] = @. 2*(t^3 - t^2) #2 is from changin from -1, 1, to 1

    #divide by 2 so each is 1 or 0 on the edges.
    S.dH[1, :] = @. (6*t^2 - 6*t) / 2
    S.dH[2, :] = @. 2*(3*t^2-4t + 1) / 2 #2 is from changin from -1, 1, to 1
    S.dH[3, :] = @. (-6t^2 + 6t) / 2
    S.dH[4, :] = @. 2*(3t^2 - 2t) / 2 #2 is from changin from -1, 1, to 1

    #these have to be divided again!
    S.ddH[1, :] = @. (12*t - 6) / 4
    S.ddH[2, :] = @. 2*(6*t-4) / 4 #2 is from changin from -1, 1, to 1
    S.ddH[3, :] = @. (-12t + 6) / 4
    S.ddH[4, :] = @. 2*(6t - 2) / 4 #2 is from changin from -1, 1, to 1

    return S

end


"""
Struct for the hermite basis functions, for 2d finite elements.
Each field has size [4, 4, rgp, θgp], dHr denotes derivative with respect to r, ddHrθ denotes ∂_r ∂_θ H.
"""
struct HB2d
    H :: Array{Float64, 4}
    dHr :: Array{Float64, 4}
    dHθ :: Array{Float64, 4}
    ddHrr :: Array{Float64, 4}
    ddHrθ :: Array{Float64, 4}
    ddHθθ :: Array{Float64, 4}
end



"""
    hermite_basis(rgp::Array{Float64}, θgp::Array{Float64})

Creates the 16 Hermite basis functions used in 2d. Creates a 1d hermite_basis for each dimension, then combines them.
"""
function hermite_basis(rgp::Array{Float64}, θgp::Array{Float64})

    #create two 1d basis'
    Sr = hermite_basis(rgp)
    Sθ = hermite_basis(θgp)
   
    #then take a tensor product to get the 2d basis
    S = combine_basis(Sr, Sθ, rgp, θgp)

    return S

end


"""
    combine_basis(Hr::Array{Float64, 2}, Hθ::Array{Float64, 2}, rgp::Array{Float64}, θgp::Array{Float64})

Combines two 1d Hermite basis' into a 2d basis.
"""
function combine_basis(Sr::HB1d, Sθ::HB1d, rgp::Array{Float64}, θgp::Array{Float64})

    #the size of each array.
    as = (4, 4, length(rgp), length(θgp))
    #define the struct to store the basis functions.
    S = HB2d(zeros(as), zeros(as), zeros(as), zeros(as), zeros(as), zeros(as))

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


    #The order of this doesn't really matter tbh.

    for i in 1:length(rgp), j in 1:length(θgp)

        #loop through the 4 nodes to consider.
        for y in 1:4, x in 1:4

            S.H[x, y, i, j] = Sr.H[x, i] * Sθ.H[y, j]
            S.dHr[x, y, i, j] = Sr.dH[x, i] * Sθ.H[y, j]
            S.dHθ[x, y, i, j] = Sr.H[x, i] * Sθ.dH[y, j]
            S.ddHrr[x, y, i, j] = Sr.ddH[x, i] * Sθ.H[y, j]
            S.ddHrθ[x, y, i, j] = Sr.dH[x, i] * Sθ.dH[y, j]
            S.ddHθθ[x, y, i, j] = Sr.H[x, i] * Sθ.ddH[y, j]
        end
    end

    return S
end



"""
Struct for the hermite basis functions, for 3d finite elements.
Each field has size [4, 4, 4, rgp, θgp, ζgp], dHr denotes derivative with respect to r, ddHrθ denotes ∂_r ∂_θ H.
"""
struct HB3d
    H :: Array{Float64, 6}
    dHr :: Array{Float64, 6}
    dHθ :: Array{Float64, 6}
    dHζ :: Array{Float64, 6}
    ddHrr :: Array{Float64, 6}
    ddHrθ :: Array{Float64, 6}
    ddHrζ :: Array{Float64, 6}
    ddHθθ :: Array{Float64, 6}
    ddHθζ :: Array{Float64, 6}
    ddHζζ :: Array{Float64, 6}
end



"""
    hermite_basis(rgp::Array{Float64}, θgp::Array{Float64}, ζgp::Array{Float64})

Creates the 64 Hermite basis functions used in 3d. Creates a 1d hermite_basis for each dimension, then combines them.
"""
function hermite_basis(rgp::Array{Float64}, θgp::Array{Float64}, ζgp::Array{Float64})

    #create 3 1d basis'
    Sr = hermite_basis(rgp)
    Sθ = hermite_basis(θgp)
    Sζ = hermite_basis(ζgp)

    #tensory product to combine 1d basis' into 3d.
    S = combine_basis(Sr, Sθ, Sζ, rgp, θgp, ζgp)
   
    return S

end


"""
    combine_basis(Hr::Array{Float64, 2}, Hθ::Array{Float64, 2}, Hζ::Array{Float64, 2}, rgp::Array{Float64}, θgp::Array{Float64}, ζgp::Array{Float64})

Combines three 1d Hermite basis' into a 3d basis.
"""
function combine_basis(Sr::HB1d, Sθ::HB1d, Sζ::HB1d, rgp::Array{Float64}, θgp::Array{Float64}, ζgp::Array{Float64})

    #size of the arrays to create.
    as = (4, 4, 4, length(rgp), length(θgp), length(ζgp))

    #creates the struct to store the functions.
    S = HB3d(zeros(as), zeros(as), zeros(as), zeros(as), zeros(as), zeros(as), zeros(as), zeros(as), zeros(as), zeros(as))

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


    for i in 1:length(rgp), j in 1:length(θgp), k in 1:length(ζgp)

        #loops over the sixteen relevant nodes, ie (0, 0, 0), (1, 0, 0) etc
        #we increment in 2's to account for derivative terms in H.
        #incrementing x,y and z will take us to the derivative at the appropriate node.
        for z in 1:4, y in 1:4, x in 1:4

            S.H[x, y, z, i, j, k] = Sr.H[x, i] * Sθ.H[y, j] * Sζ.H[z, k]
            S.dHr[x, y, z, i, j, k] = Sr.dH[x, i] * Sθ.H[y, j] * Sζ.H[z, k]
            S.dHθ[x, y, z, i, j, k] = Sr.H[x, i] * Sθ.dH[y, j] * Sζ.H[z, k]
            S.dHζ[x, y, z, i, j, k] = Sr.H[x, i] * Sθ.H[y, j] * Sζ.dH[z, k]
            S.ddHrr[x, y, z, i, j, k] = Sr.ddH[x, i] * Sθ.H[y, j] * Sζ.H[z, k]
            S.ddHrθ[x, y, z, i, j, k] = Sr.dH[x, i] * Sθ.dH[y, j] * Sζ.H[z, k]
            S.ddHrζ[x, y, z, i, j, k] = Sr.dH[x, i] * Sθ.H[y, j] * Sζ.dH[z, k]
            S.ddHθθ[x, y, z, i, j, k] = Sr.H[x, i] * Sθ.ddH[y, j] * Sζ.H[z, k]
            S.ddHθζ[x, y, z, i, j, k] = Sr.H[x, i] * Sθ.dH[y, j] * Sζ.dH[z, k]
            S.ddHζζ[x, y, z, i, j, k] = Sr.H[x, i] * Sθ.H[y, j] * Sζ.ddH[z, k]

        end
    end
    
    return S
end


