
"""
Struct storing the coeffiecients of the qfm surfaces.
Based on expansion 
ψ(φ) = ∑ ψ_n^c cos(φ) + ψ_n^s sin(φ)
θ(φ) = θ_0^c + bφ/a + ∑ θ_n^c cos(φ) + θ_n^s sin(φ)
"""
struct QFMSurfaceT
    rational :: Tuple{Int64, Int64} #(a, b)
    s :: Float64 #surface label
    ψcos :: Array{Float64, 2}
    θsin :: Array{Float64, 2}
    ψsin :: Array{Float64, 2}
    θcos :: Array{Float64, 2}
end


"""
Struct for storing the temporary values defining each surface used in coord_transform!
"""
struct TempSurfT
    ψcos :: Array{Float64, 2} 
    θsin :: Array{Float64, 2}
    dψcosds :: Array{Float64, 2}
    dθsinds :: Array{Float64, 2}
    d2ψcosdsds :: Array{Float64, 2}
    d2θsindsds :: Array{Float64, 2}
    α :: Array{Float64, 2}
    cosα :: Array{Float64, 2}
    sinα :: Array{Float64, 2}
    function TempSurfT(M::Int64, N::Int64) 
        dim1 = M + 1
        dim2 = 2*N + 1
        new(zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2))
    end
end


"""
Struct storing the interpolations used to create surfaces in between the qfm surfaces found.
Derivatives are stored individually for extrapolation to work.
"""
struct SurfaceITPT
    M :: Int64
    N :: Int64
    ψcos_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    θsin_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    dψcos_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    dθsin_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    d2ψcos_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    d2θsin_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
end



