

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
    hermite_basis(rgp::Array{Float64}, θgp::Array{Float64}, ζgp::Array{Float64})

Creates the 64 Hermite basis functions used in 3d. Creates a 1d hermite_basis for each dimension, then combines them.
"""
function hermite_basis(rgp::Array{Float64}, θgp::Array{Float64}, ζgp::Array{Float64})

    Hr, dHr, ddHr = hermite_basis(rgp)
    Hθ, dHθ, ddHθ = hermite_basis(θgp)
    Hζ, dHζ, ddHζ = hermite_basis(ζgp)

    #probably could be more efficient, but this is only done once so shouldn't matter.
    S = combine_basis(Hr, Hθ, Hζ, rgp, θgp, ζgp)
    dSr = combine_basis(dHr, Hθ, Hζ, rgp, θgp, ζgp)
    dSθ = combine_basis(Hr, dHθ, Hζ, rgp, θgp, ζgp)
    dSζ = combine_basis(Hr, Hθ, dHζ, rgp, θgp, ζgp)
    ddSrr = combine_basis(ddHr, Hθ, Hζ, rgp, θgp, ζgp)
    ddSrθ = combine_basis(dHr, dHθ, Hζ, rgp, θgp, ζgp)
    ddSrζ = combine_basis(dHr, Hθ, dHζ, rgp, θgp, ζgp)
    ddSθθ = combine_basis(Hr, ddHθ, Hζ, rgp, θgp, ζgp)
    ddSθζ = combine_basis(Hr, dHθ, dHζ, rgp, θgp, ζgp)
    ddSζζ = combine_basis(Hr, Hθ, ddHζ, rgp, θgp, ζgp)
    
    return S, dSr, dSθ, dSζ, ddSrr, ddSrθ, ddSrζ, ddSθθ, ddSθζ, ddSζζ

end

"""
    combine_basis(Hr::Array{Float64, 2}, Hθ::Array{Float64, 2}, Hζ::Array{Float64, 2}, rgp::Array{Float64}, θgp::Array{Float64}, ζgp::Array{Float64})

Combines three 1d Hermite basis' into a 3d basis.
"""
function combine_basis(Hr::Array{Float64, 2}, Hθ::Array{Float64, 2}, Hζ::Array{Float64, 2}, rgp::Array{Float64}, θgp::Array{Float64}, ζgp::Array{Float64})


    S = zeros(Float64, 4, 4, 4, length(rgp), length(θgp), length(ζgp))

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

            S[x, y, z, i, j, k] = Hr[x, i] * Hθ[y, j] * Hζ[z, k]
        end
    end
    
    return S
end




"""
    hermite_basis(rgp::Array{Float64}, θgp::Array{Float64})

Creates the 16 Hermite basis functions used in 2d. Creates a 1d hermite_basis for each dimension, then combines them.
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


"""
    combine_basis(Hr::Array{Float64, 2}, Hθ::Array{Float64, 2}, rgp::Array{Float64}, θgp::Array{Float64})

Combines two 1d Hermite basis' into a 2d basis.
"""
function combine_basis(Hr::Array{Float64, 2}, Hθ::Array{Float64, 2}, rgp::Array{Float64}, θgp::Array{Float64})


    S = zeros(4, length(rgp), 4, length(θgp))

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

            S[x, y, i, j] = Hr[x, i] * Hθ[y, j]
        end
    end

    return S
end


