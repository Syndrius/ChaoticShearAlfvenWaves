
"""
Struct storing the grids.
# Fields
- rd::RDataT Struct storing the finite elements radial grid information.
- pmd::ModeDadaT Struct storing the spectral method poloidal grid information.
- tmd::ModeDadaT Struct storing the spectral method toroidal grid information.
"""
@kwdef struct GridsT
    rd :: RDataT
    pmd :: ModeDataT
    tmd :: ModeDataT
end


"""
Constructor for struct storing the grids.
# Fields
- rgrid::ArrayFloat64 Array of radial points.
- mstart::Int64 First poloidal mode.
- mcount::Int64 Number of poloidal modes.
- mincr::Int64=1 Gap between poloidal modes.
- nstart::Int64 First toroidal mode.
- ncount::Int64 Number of toroidal modes.
- nincr::Int64=1 Gap between toroidal modes.
- g_quad::Int64=4 Gaussian quadrature for finite elements integration.
- f_quad::Int64=3 Fourier quadrature for defining angular grids.
"""
function init_grids(; rgrid::Array{Float64}, mstart::Int64, mcount::Int64, mincr::Int64=1, nstart::Int64, ncount::Int64, nincr::Int64=1, g_quad::Int64=4, f_quad::Int64=3)

    rd = RDataT(collect(rgrid), g_quad, length(rgrid))
    pmd = ModeDataT(mstart, mcount, mincr, f_quad)
    tmd = ModeDataT(nstart, ncount, nincr, f_quad)

    return GridsT(rd=rd, pmd=pmd, tmd=tmd)

end