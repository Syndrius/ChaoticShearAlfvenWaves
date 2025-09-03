
"""
Struct storing the coeffiecients of the qfm surfaces.
Based on expansion 
r(ζ) = ∑ r_n^c cos(ζ) + r_n^s sin(ζ)
θ(ζ) = θ_0^c + bζ/a + ∑ θ_n^c cos(ζ) + θ_n^s sin(ζ)
"""
struct QFMSurfaceT
    rational :: Tuple{Int64, Int64} #(a, b)
    s :: Float64 #surface label
    rcos :: Array{Float64, 2}
    θsin :: Array{Float64, 2}
    rsin :: Array{Float64, 2}
    θcos :: Array{Float64, 2}
end

#name subject to change.
"""
Struct for storing the temporary values defining each surface used in coord_transform!
"""
struct TempSurfT
    rcos :: Array{Float64, 2} 
    θsin :: Array{Float64, 2}
    drcosds :: Array{Float64, 2}
    dθsinds :: Array{Float64, 2}
    d2rcosdsds :: Array{Float64, 2}
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
    rcos_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    θsin_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    drcos_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    dθsin_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    d2rcos_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    d2θsin_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
end



