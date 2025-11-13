#name and function of this file is subject to change!
#with this layout of grids, we could pass in a single ξ tuple object

#maybe ξ should be a static array!
function create_basis(grids::FSSGridsT, ξ::Array{Float64})

    return create_hermite_basis(ξ)
end

function create_basis(grids::FFSGridsT, ξ1::Array{Float64}, ξ2::Array{Float64})

    S1 = create_hermite_basis(ξ1)
    S2 = create_hermite_basis(ξ2)

    return combine_basis(S1, S2, ξ1, ξ2) #needs to change
end

function create_basis(grids::FFFGridsT, ξ1::Array{Float64}, ξ2::Array{Float64}, ξ3::Array{Float64})

    S1 = create_hermite_basis(ξ1)
    S2 = create_hermite_basis(ξ2)
    S3 = create_hermite_basis(ξ3)

    #don't really need the ξ's...
    return combine_basis(S1, S2, S3, ξ1, ξ2, ξ3)
end

function create_basis(grids::FSSGridsT, ξ::Array{Float64}, mlist, nlist)

    S1 = create_hermite_basis(ξ)
    S2 = create_spectral_basis(mlist)
    S3 = create_spectral_basis(nlist)

    #even this is a waste, the info is within S right?
    return combine_basis(S1, S2, S3, length(ξ), length(mlist), length(nlist))
end
