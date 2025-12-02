
"""
    function compute_P!(P::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, B::BFieldT, met::MetT, n::Float64, tm::TM)

Computes the P matrix for the weak form at a single coordinate.
"""
function compute_P!(P::SubArray{ComplexF64, 2, Array{ComplexF64, 5}}, B::BFieldT, met::MetT, n::Float64, tm::TM)

   
    #compute the laplacian like term
    Tl!(P, B, met, tm.C, tm.D)
    
    #compute the current term.
    Tj!(P, B, met, tm.Γ, tm.dΓ, tm.K)
    

end

