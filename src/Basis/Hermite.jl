
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
Each field has size [4, 4, x1gp, x2gp], dHx1 denotes derivative with respect to x1, ddHx1x2 denotes ∂_x1 ∂_x2 H.
"""
struct HB2d
    H :: Array{Float64, 4}
    dHx1 :: Array{Float64, 4}
    dHx2 :: Array{Float64, 4}
    ddHx1x1 :: Array{Float64, 4}
    ddHx1x2 :: Array{Float64, 4}
    ddHx2x2 :: Array{Float64, 4}
end



"""
    hermite_basis(x1gp::Array{Float64}, x2gp::Array{Float64})

Creates the 16 Hermite basis functions used in 2d. Creates a 1d hermite_basis for each dimension, then combines them.
"""
function hermite_basis(x1gp::Array{Float64}, x2gp::Array{Float64})

    #create two 1d basis'
    Sx1 = hermite_basis(x1gp)
    Sx2 = hermite_basis(x2gp)
   
    #then take a tensor product to get the 2d basis
    S = combine_basis(Sx1, Sx2, x1gp, x2gp)

    return S

end


"""
    combine_basis(Hx1::Array{Float64, 2}, Hx2::Array{Float64, 2}, x1gp::Array{Float64}, x2gp::Array{Float64})

Combines two 1d Hermite basis' into a 2d basis.
"""
function combine_basis(Sx1::HB1d, Sx2::HB1d, x1gp::Array{Float64}, x2gp::Array{Float64})

    #the size of each array.
    as = (4, 4, length(x1gp), length(x2gp))
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

    for i in 1:length(x1gp), j in 1:length(x2gp)

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
struct HB3d
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
    hermite_basis(x1gp::Array{Float64}, x2gp::Array{Float64}, x3gp::Array{Float64})

Creates the 64 Hermite basis functions used in 3d. Creates a 1d hermite_basis for each dimension, then combines them.
"""
function hermite_basis(x1gp::Array{Float64}, x2gp::Array{Float64}, x3gp::Array{Float64})

    #create 3 1d basis'
    Sx1 = hermite_basis(x1gp)
    Sx2 = hermite_basis(x2gp)
    Sx3 = hermite_basis(x3gp)

    #tensory product to combine 1d basis' into 3d.
    S = combine_basis(Sx1, Sx2, Sx3, x1gp, x2gp, x3gp)
   
    return S

end


"""
    combine_basis(Hx1::Array{Float64, 2}, Hx2::Array{Float64, 2}, Hx3::Array{Float64, 2}, x1gp::Array{Float64}, x2gp::Array{Float64}, x3gp::Array{Float64})

Combines three 1d Hermite basis' into a 3d basis.
"""
function combine_basis(Sx1::HB1d, Sx2::HB1d, Sx3::HB1d, x1gp::Array{Float64}, x2gp::Array{Float64}, x3gp::Array{Float64})

    #size of the arrays to create.
    as = (4, 4, 4, length(x1gp), length(x2gp), length(x3gp))

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


    for i in 1:length(x1gp), j in 1:length(x2gp), k in 1:length(x3gp)

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


