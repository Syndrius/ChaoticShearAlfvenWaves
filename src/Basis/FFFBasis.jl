
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
    create_local_basis!(Φ::Array{ComplexF64, 7}, S::Array{Float64, 6}, dSr::Array{Float64, 6}, dSθ::Array{Float64, 6}, dSζ::Array{Float64, 6}, ddSrr::Array{Float64, 6}, ddSrθ::Array{Float64, 6}, ddSrζ::Array{Float64, 6}, ddSθθ::Array{Float64, 6}, ddSθζ::Array{Float64, 6}, ddSζζ::Array{Float64, 6}, m::Int64, n::Int64, dr::Float64, dθ::Float64, dζ::Float64)

Modifies the matrix storing the basis functions, Φ, to reflect the local derivatives being taken. Includes a phase factor modification for θ and ζ.
"""
function create_local_basis!(Φ::Array{ComplexF64, 7}, S::Array{Float64, 6}, dSr::Array{Float64, 6}, dSθ::Array{Float64, 6}, dSζ::Array{Float64, 6}, ddSrr::Array{Float64, 6}, ddSrθ::Array{Float64, 6}, ddSrζ::Array{Float64, 6}, ddSθθ::Array{Float64, 6}, ddSθζ::Array{Float64, 6}, ddSζζ::Array{Float64, 6}, m::Int64, n::Int64, dr::Float64, dθ::Float64, dζ::Float64)

    #TODO Change the inputs, they are stupid af, Should collect all S into one!

    rjac = dr / 2
    θjac = dθ / 2
    ζjac = dζ / 2
    
    #Φ_s
    @. Φ[:, :, :, 1, :, :, :] = dSr / rjac 
    #Φ_θ
    @. Φ[:, :, :, 2, :, :, :] = dSθ / θjac + S * 1im * m
    #Φ_ζ
    @. Φ[:, :, :, 3, :, :, :] = dSζ / ζjac + S * 1im * n
    #Φ_ss
    @. Φ[:, :, :, 4, :, :, :] = ddSrr / rjac^2 
    #Φ_sθ
    @. Φ[:, :, :, 5, :, :, :] = ddSrθ / (rjac * θjac) + dSr * 1im * m / rjac
    #Φ_sζ   
    @. Φ[:, :, :, 6, :, :, :] = ddSrζ / (rjac * ζjac) + dSr * 1im * n / rjac 
    #Φ_θθ
    @. Φ[:, :, :, 7, :, :, :] = ddSθθ / θjac^2 + 2 * dSθ * 1im * m / θjac - m^2 * S
    #Φ_θζ
    @. Φ[:, :, :, 8, :, :, :] = ddSθζ / (θjac * ζjac) + dSζ * 1im * m / ζjac + dSθ * 1im * n / θjac - m * n * S
    #Φ_ζζ
    @. Φ[:, :, :, 9, :, :, :] = ddSζζ / ζjac^2 + 2 * dSζ * 1im * n / ζjac - n^2 * S

end

        