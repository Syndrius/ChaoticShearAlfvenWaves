
#define the order of the shape functions.
const grid_id = [0, 0, 1, 1]
const basis_id = [0, 1, 0, 1]


"""
Struct for the hermite basis functions, for 1d finite elements.

### Fields
- H::Array{Float64, 2} - Shape functions, size [4, gp]
- dH::Array{Float64, 2} - Derivatives, size [4, gp]
- ddH::Array{Float64, 2} - Second derivatives, size [4, gp]
"""
struct HB1dT
    H :: Array{Float64, 2}
    dH :: Array{Float64, 2}
    ddH :: Array{Float64, 2}
end


"""
    hermite_basis(gp::Array{Float64})

Creates the four Hermite basis functions, taken from wikipedia. Returns a 4xgp matrix for the 0th, 1st and 2nd derivative.
"""
function hermite_basis(ξ::Array{Float64})

    gp = length(ξ) #number of grid points in the gaussian quadrature

    #creates a struct to store the hermiete shape functions, their derivative and second derivative.
    S = HB1dT(zeros(4, gp), zeros(4, gp), zeros(4, gp))

    ##############
    # node structure per element
    #element defined so locally x = [-1, 1]
    # -1 -------- 1

    #each of the basis functions has either a value of 1 at one edge or a derivative of 1 at one edge
    #and then value of zero at the other three options.

    #ie first basis functions has H(-1) = 1, H'(-1) = 0, H(1) = 0, H'(1) = 0

    #node 1 H(-1) = 1
    S.H[1, :] = h00.(ξ)

    #node 1 H'(-1) = 1
    S.H[2, :] = h10.(ξ)

    #node 2 H(1) = 1
    S.H[3, :] = h01.(ξ)

    #node 2 H'(1) = 1
    S.H[4, :] = h11.(ξ)

    S.dH[1, :] = dh00.(ξ)
    S.dH[2, :] = dh10.(ξ)
    S.dH[3, :] = dh01.(ξ)
    S.dH[4, :] = dh11.(ξ)

    S.ddH[1, :] = ddh00.(ξ)
    S.ddH[2, :] = ddh10.(ξ)
    S.ddH[3, :] = ddh01.(ξ)
    S.ddH[4, :] = ddh11.(ξ)

    return S

end


"""
Struct for the hermite basis functions, for 2d finite elements.
Each field has size [4, 4, x1gp, x2gp], dHx1 denotes derivative with respect to x1, ddHx1x2 denotes ∂_x1 ∂_x2 H.
"""
struct HB2dT
    H :: Array{Float64, 4}
    dHx1 :: Array{Float64, 4}
    dHx2 :: Array{Float64, 4}
    ddHx1x1 :: Array{Float64, 4}
    ddHx1x2 :: Array{Float64, 4}
    ddHx2x2 :: Array{Float64, 4}
end



"""
    hermite_basis(ξ1::Array{Float64}, ξ2::Array{Float64})

Creates the 16 Hermite basis functions used in 2d. Creates a 1d hermite_basis for each dimension, then combines them.
"""
function hermite_basis(ξ1::Array{Float64}, ξ2::Array{Float64})

    #create two 1d basis'
    Sx1 = hermite_basis(ξ1)
    Sx2 = hermite_basis(ξ2)
   
    #then take a tensor product to get the 2d basis
    S = combine_basis(Sx1, Sx2, ξ1, ξ2)

    return S

end


"""
    combine_basis(Hx1::Array{Float64, 2}, Hx2::Array{Float64, 2}, ξ1::Array{Float64}, ξ2::Array{Float64})

Combines two 1d Hermite basis' into a 2d basis.
"""
function combine_basis(Sx1::HB1dT, Sx2::HB1dT, ξ1::Array{Float64}, ξ2::Array{Float64})

    #the size of each array.
    as = (4, 4, length(ξ1), length(ξ2))
    #define the struct to store the basis functions.
    S = HB2dT(zeros(as), zeros(as), zeros(as), zeros(as), zeros(as), zeros(as))

    #recall from 1d
    ##############
    # node structure per element

    # -1 -------- 1

    #H[1] -> H(-1) = 1
    #H[2] -> H'(-1) = 1
    #H[3] -> H(1) = 1
    #H[4] -> H'(1) = 1


    ############### 
    #node labels for each element!

    # (-1, 1) -------- (1, 1)
    #  |               |
    #  |               |
    #  |               |
    # (-1, -1) -------- (1, -1)

    for i in 1:length(ξ1), j in 1:length(ξ2)

        #loop through the 4 nodes to consider.
        for y in 1:4, x in 1:4

            S.H[x, y, i, j] = Sx1.H[x, i] * Sx2.H[y, j]
            S.dHx1[x, y, i, j] = Sx1.dH[x, i] * Sx2.H[y, j]
            S.dHx2[x, y, i, j] = Sx1.H[x, i] * Sx2.dH[y, j]
            S.ddHx1x1[x, y, i, j] = Sx1.ddH[x, i] * Sx2.H[y, j]
            S.ddHx1x2[x, y, i, j] = Sx1.dH[x, i] * Sx2.dH[y, j]
            S.ddHx2x2[x, y, i, j] = Sx1.H[x, i] * Sx2.ddH[y, j]
        end
    end

    return S
end



"""
Struct for the hermite basis functions, for 3d finite elements.
Each field has size [4, 4, 4, x1gp, x2gp, x3gp], dHx1 denotes derivative with respect to x1, ddHx1x2 denotes ∂_x1 ∂_x2 H.
"""
struct HB3dT
    H :: Array{Float64, 6}
    dHx1 :: Array{Float64, 6}
    dHx2 :: Array{Float64, 6}
    dHx3 :: Array{Float64, 6}
    ddHx1x1 :: Array{Float64, 6}
    ddHx1x2 :: Array{Float64, 6}
    ddHx1x3 :: Array{Float64, 6}
    ddHx2x2 :: Array{Float64, 6}
    ddHx2x3 :: Array{Float64, 6}
    ddHx3x3 :: Array{Float64, 6}
end



"""
    hermite_basis(ξ1::Array{Float64}, x2gp::Array{Float64}, ξ3::Array{Float64})

Creates the 64 Hermite basis functions used in 3d. Creates a 1d hermite_basis for each dimension, then combines them.
"""
function hermite_basis(ξ1::Array{Float64}, ξ2::Array{Float64}, ξ3::Array{Float64})

    #create 3 1d basis'
    Sx1 = hermite_basis(ξ1)
    Sx2 = hermite_basis(ξ2)
    Sx3 = hermite_basis(ξ3)

    #tensory product to combine 1d basis' into 3d.
    S = combine_basis(Sx1, Sx2, Sx3, ξ1, ξ2, ξ3)
   
    return S

end


"""
    combine_basis(Hx1::Array{Float64, 2}, Hx2::Array{Float64, 2}, Hx3::Array{Float64, 2}, ξ1::Array{Float64}, ξ2::Array{Float64}, ξ3::Array{Float64})

Combines three 1d Hermite basis' into a 3d basis.
"""
function combine_basis(Sx1::HB1dT, Sx2::HB1dT, Sx3::HB1dT, ξ1::Array{Float64}, ξ2::Array{Float64}, ξ3::Array{Float64})

    #size of the arrays to create.
    as = (4, 4, 4, length(ξ1), length(ξ2), length(ξ3))

    #creates the struct to store the functions.
    S = HB3dT(zeros(as), zeros(as), zeros(as), zeros(as), zeros(as), zeros(as), zeros(as), zeros(as), zeros(as), zeros(as))

    #recall from 1d
    ##############
    # node structure per element

    # -1 -------- 1

    #H[1] -> H(-1) = 1
    #H[2] -> H'(-1) = 1
    #H[3] -> H(1) = 1
    #H[4] -> H'(1) = 1


    ############### 
    #node labels for each element!

    # (-1, 1) -------- (1, 1)
    #  |               |
    #  |               |
    #  |               |
    # (-1, -1) -------- (1, -1)

    for i in 1:length(ξ1), j in 1:length(ξ2), k in 1:length(ξ3)

        #loops over the sixteen relevant nodes, ie (0, 0, 0), (1, 0, 0) etc
        #we increment in 2's to account for derivative terms in H.
        #incrementing x,y and z will take us to the derivative at the appropriate node.
        for z in 1:4, y in 1:4, x in 1:4

            S.H[x, y, z, i, j, k] = Sx1.H[x, i] * Sx2.H[y, j] * Sx3.H[z, k]
            S.dHx1[x, y, z, i, j, k] = Sx1.dH[x, i] * Sx2.H[y, j] * Sx3.H[z, k]
            S.dHx2[x, y, z, i, j, k] = Sx1.H[x, i] * Sx2.dH[y, j] * Sx3.H[z, k]
            S.dHx3[x, y, z, i, j, k] = Sx1.H[x, i] * Sx2.H[y, j] * Sx3.dH[z, k]
            S.ddHx1x1[x, y, z, i, j, k] = Sx1.ddH[x, i] * Sx2.H[y, j] * Sx3.H[z, k]
            S.ddHx1x2[x, y, z, i, j, k] = Sx1.dH[x, i] * Sx2.dH[y, j] * Sx3.H[z, k]
            S.ddHx1x3[x, y, z, i, j, k] = Sx1.dH[x, i] * Sx2.H[y, j] * Sx3.dH[z, k]
            S.ddHx2x2[x, y, z, i, j, k] = Sx1.H[x, i] * Sx2.ddH[y, j] * Sx3.H[z, k]
            S.ddHx2x3[x, y, z, i, j, k] = Sx1.H[x, i] * Sx2.dH[y, j] * Sx3.dH[z, k]
            S.ddHx3x3[x, y, z, i, j, k] = Sx1.H[x, i] * Sx2.H[y, j] * Sx3.ddH[z, k]

        end
    end
    
    return S
end

