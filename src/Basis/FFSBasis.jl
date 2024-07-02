


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



"""
    create_local_basis!(Φ::Array{ComplexF64, 5}, S::Array{Float64, 4}, dSr::Array{Float64, 4}, dSθ::Array{Float64, 4}, ddSrr::Array{Float64, 4}, ddSrθ::Array{Float64, 4}, ddSθθ::Array{Float64, 4}, m::Int64, n::Int64, dr::Float64, dθ::Float64)

Modifies the matrix storing the basis functions, Φ, to reflect the local derivatives being taken. Includes a phase factor modification for θ.
"""
function create_local_basis!(Φ::Array{ComplexF64, 5}, S::Array{Float64, 4}, dSr::Array{Float64, 4}, dSθ::Array{Float64, 4}, ddSrr::Array{Float64, 4}, ddSrθ::Array{Float64, 4}, ddSθθ::Array{Float64, 4}, m::Int64, n::Int64, dr::Float64, dθ::Float64)

    rjac = dr / 2
    θjac = dθ / 2
    
    #Φ_s
    @. Φ[:, :, 1, :, :] = dSr / rjac 
    #Φ_θ
    @. Φ[:, :, 2, :, :] = dSθ / θjac + S * 1im * m
    #Φ_ζ
    @. Φ[:, :, 3, :, :] = S * n * 1im 
    #Φ_ss
    @. Φ[:, :, 4, :, :] = ddSrr / rjac^2 
    #Φ_sθ
    @. Φ[:, :, 5, :, :] = ddSrθ / (rjac * θjac) + dSr * 1im * m / rjac
    #Φ_sζ   
    @. Φ[:, :, 6, :, :] = dSr * n * 1im / rjac 
    #Φ_θθ
    @. Φ[:, :, 7, :, :] = ddSθθ / θjac^2 + 2 * dSθ * 1im * m / θjac - m^2 * S
    #Φ_θζ
    @. Φ[:, :, 8, :, :] = 1im * n * (dSθ / θjac + S * 1im * m)
    #Φ_ζζ
    @. Φ[:, :, 9, :, :] = S * (-n^2) 

end



