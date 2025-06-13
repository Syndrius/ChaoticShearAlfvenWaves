

"""
    f2r(ψ)

Converts the flux to radius.
"""
function f2r(ψ)
    #assumes B0=1
    return sqrt(2*ψ)
end


"""
    r2f(r)

Converts radius to toroidal flux.
"""
function r2f(r)
    #assumes B0=1
    return r^2/2
end


function convert_island(isl::RadIslandT)

    #need to make sure this actually covers everything first
    #this also will only work for a fully instantiated island.
    qp = isl.qp / isl.r0
    A = (isl.w/4)^2 * qp / isl.q0^2
    ψ0 = isl.r0^2 / 2
    return PsiIslandT(isl.m0, isl.n0, A, isl.q0, qp, ψ0)

end


function convert_island(isl::FluxIslandT)
    #TODO

end
