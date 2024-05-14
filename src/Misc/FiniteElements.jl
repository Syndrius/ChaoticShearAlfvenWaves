
"""
Struct storing the data on the radial grid for finite elements.

# Fields
- grid::Array{Float64} Array storing the grid.
- gp::Int64=4 Number of Gaussian quadrature points to use for numerical integration, defaults to 4.
- N::Int64 Size of the grid.
"""
@kwdef struct RDataT
    grid :: Array{Float64}
    gp :: Int64 = 4
    N :: Int64
end

"""
Two arrays that distinguish the different Hermite basis functions, used for converting between the grid and the appropriate points indices in the matrices.
"""
const grid_id = [0, 0, 1, 1]
const basis_id = [0, 1, 0, 1]


"""
Creates the four Hermite basis functions, taken from wikipedia. Returns a 4xgp matrix for the 0th, 1st and 2nd derivative.

# Args
- gp::Array{Float64} Array of the Gaussian quadrature points.
"""
function hermite_basis(gp::Array{Float64})
    t = @. (gp + 1)/(2) #converts to correct range for spline
    H = zeros(4, length(gp))
    dH = zeros(4, length(gp))
    ddH = zeros(4, length(gp))

    H[1, :] = @. 2*t^3 - 3*t^2 + 1
    H[2, :] = @. 2*(t^3-2t^2 + t) #2 is from changin from -1, 1, to 1
    H[3, :] = @. -2t^3 + 3t^2
    H[4, :] = @. 2*(t^3 - t^2) #2 is from changin from -1, 1, to 1
    #display(H)

    #these all seem to need to be divided by 2, not sure why
    #think there is an extra condition on these being 1/0 on the edges.
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
Converts the local grid where the finite elements are defined to the global coordinates.

# Args
- node::Int64 Current node in the finite elements grid.
- ξ::Array{Float64} Local grid, -1≤ξ≤1.
- grid::Array{Float64} Global grid.
"""
function local_to_global(node::Int64, ξ::Array{Float64}, grid::Array{Float64})

    dr = grid[node+1] - grid[node]

    #first we map the (-1, 1) local coords to (0, 1)
    mp = @. (ξ+1) / 2

    #then we sale by the width of grid points
    mp = mp .* dr

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    rglobal = mp .+ grid[node]

    
    return rglobal, dr
end 






"""
Computes the numerical integration required for finite elements via ∫f(x)dx ≈ ∑w_i f(x_i).
#this function contracts the matrix with the basis we have made and integrates using gauss quadrature.
#always 9x9 now.
#these are extra stupid types as we are taking views, we have just copied the error message, this may be a bad idea!

# Args These are stupid af and probably subject to change.
"""
function gauss_integrate(test_vec::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, trial_vec::SubArray{ComplexF64, 2, Array{ComplexF64, 3}}, mat::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, wg::Array{Float64}, jac::Float64, ngp::Int64)#::ComplexF64

    #this function is taking up a lot of time when nm, nn get large, but I think that is more to do with the size of the loops that call this function tbh.

    res = 0.0 + 0.0im
    for k in 1:ngp
    
        scale = wg[k] * jac
        for j in 1:9

            for i in 1:9
                #significantly faster to have * jac here not later!
                res += @inbounds test_vec[i, k] * mat[i, j, k] * trial_vec[j, k] * scale
            #display(res)
            end
        end
    end
    #res *= jac
    return res #* dr / 2 
end


#=
#not used anymore, checked if it was faster, but it is not, kinda surprised tbh!
function combined_integrate(resI, resW, test_vec, trial_vec, I, W, wg, jac, ngp, left_dim, right_dim)#::Tuple{ComplexF64, ComplexF64}
    #this one is a wee bit slower for some reason, has more gc time..
    resI = 0.0 + 0.0im
    resW = 0.0 + 0.0im
    for i in 1:left_dim

        for j in 1:right_dim

            for k in 1:ngp

                resI += @inbounds @views test_vec[i, k] * I[i, j, k] * trial_vec[j, k] * wg[k] * jac
                resW += @inbounds @views test_vec[i, k] * W[i, j, k] * trial_vec[j, k] * wg[k] * jac
            #display(res)
            end
        end
    end
    #return (resI * dr / 2, resW * dr / 2)

end
=#