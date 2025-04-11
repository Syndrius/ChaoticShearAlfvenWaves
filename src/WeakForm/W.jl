
"""
    function compute_W!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, B::BFieldT, met::MetT, n::Float64, ωcap2::Float64, tm::TM)

Computes the W matrix for the weak form at a single coordinate.
"""
function compute_W!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, B::BFieldT, met::MetT, n::Float64, ωcap2::Float64, tm::TM)


   
    #compute the laplacian like term
    Tl!(W, B, met, tm.C, tm.D, tm.T)
    
    #compute the current term.
    Tj!(W, B, met, tm.Γ, tm.dΓ, tm.K)
    

    #TODO
    #work in progress. This is the BAE like term, to make (0, 0) harmonic more physical.
    #W[1:3, 1:3] += tm.D .* (ωcap2 * n *  met.J / B.mag_B^2)

end


"""
    function compute_isl_W!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, B::BFieldT, met::MetT, tm::TM)

Computes the W matrix in the case of island coordinates. In this case, the current term is excluded.
"""
function compute_isl_W!(W::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, B::BFieldT, met::MetT, tm::TM)

    #only the laplace term is used in this case.
    Tl!(W, B, met, tm.C, tm.D, tm.T)
    
end

