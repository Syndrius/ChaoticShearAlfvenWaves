

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