
#TODO Split this file up
#this might mean the end of Misc, as we can probably have an integration, basis and indexing module now.
#may need to move some of the gauss quadrature stuff to integration.




########################################
#Legacy stuff for combining the basis functions.


#no longer used!
const grid_id_r = [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1]
const grid_id_θ = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
const basis_id_r = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
const basis_id_θ = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]



#creates the local basis for 2d fem case.
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



#combines two 1d Hermite basis' into a 2d basis.
function og_combine_basis(Hr, Hθ, rgp, θgp)


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




function combine_basis_comb(Hr, Hθ, rgp, θgp)


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

    #we need to do this in a smarter way!

    for i in 1:length(rgp), j in 1:length(θgp)

        #loop through the 4 nodes to consider.
        for y in 1:4, x in 1:4

            #const grid_id = [0, 0, 1, 1]
            #const basis_id = [0, 1, 0, 1]

            ind = 1 + basis_id[x] + 2*basis_id[y] + 4*grid_id[x] + 8*grid_id[y]
            S[ind, i, j] = Hr[x, i] * Hθ[y, j]
        end
    end

    return S
end

