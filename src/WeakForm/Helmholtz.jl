

#most of this is obviously not being used, but this keeps the form the same so construct is not modified.
function W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, B::BFieldT, met::MetT, prob::Helmholtz.TestProblemT, r::Array{Float64}, θ::AbstractArray, ζ::AbstractArray, tm::TM)

    for k in 1:length(ζ), j in 1:length(θ), i in 1:length(r)


        #cylindrical case.
        #may want to consider swapping I and W
        #and solving for 1/k^2
        I[3, 3, i, j, k] = 1.0

        W[6, 6, i, j, k] = 1.0

        W[3, 6, i, j, k] = -1.0 / r[i]

        W[8, 8, i, j, k] = 1.0 / r[i]^2

        W[9, 9, i, j, k] = 1.0


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
