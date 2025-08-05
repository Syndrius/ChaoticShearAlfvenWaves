

#most of this is obviously not being used, but this keeps the form the same so construct is not modified.
function W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, B::BFieldT, met::MetT, prob::Helmholtz.TestProblemT, r::Array{Float64}, θ::AbstractArray, ζ::AbstractArray, tm::TM)

    for k in 1:length(ζ), j in 1:length(θ), i in 1:length(r)


        #cylindrical case.
        #artifically adds extra second derivative to ζ, 
        #doesn't really work! Think periodic boundary conditions cook this a wee bit.
        #I[3, 3, i, j, k] = 1.0

        #W[6, 6, i, j, k] = 1.0

        #W[3, 6, i, j, k] = -1.0 / r[i]

        #W[8, 8, i, j, k] = 1.0 / r[i]^2

        #W[9, 9, i, j, k] = 1.0

        #assuming the extra zeroth derivative is added to 10th index
        I[10, 10, i, j, k] = 1.0

        W[1, 1, i, j, k] = 1.0

        W[10, 1, i, j, k] = -1.0 / r[i]

        W[2, 2, i, j, k] = 1.0 / r[i]^2

        W[3, 3, i, j, k] = 1.0


        #cartesian case, requires serious modification to make work.
        #ζ
        #I[3, 3, i, j, k] = 1.0
        #rζ
        #W[6, 6, i, j, k] = 1.0
        #θζ
        #W[8, 8, i, j, k] = 1.0
        #ζζ
        #W[9, 9, i, j, k] = 1.0
    end


end
