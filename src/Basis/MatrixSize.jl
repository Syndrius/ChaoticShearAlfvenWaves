"""
    matrix_size(grids::FSSGridsT)

Computes the dimension of the matrix.
"""
function matrix_size(grids::ContGridsT)

    return grids.θ.N * grids.ζ.N
end


"""
    matrix_size(grids::FSSGridsT)

Computes the dimension of the matrix.
"""
function matrix_size(grids::FSSGridsT)

    return 2 * grids.r.N * grids.θ.N * grids.ζ.N
end


"""
    matrix_size(grids::FFSGridsT)

Computes the dimension of the matrix.
"""
function matrix_size(grids::FFSGridsT)

    return 4 * grids.r.N * grids.θ.N * grids.ζ.N
end


"""
    matrix_size(grids::FFFGridsT)

Computes the dimension of the matrix.
"""
function matrix_size(grids::FFFGridsT)

    return 8 * grids.r.N * grids.θ.N * grids.ζ.N
end

"""
    local_matrix_size(grids::FSSGridsT)

Computes the size of the local I, W matrices needed at each grid point.
"""
function local_matrix_size(grids::ContGridsT)

    return zeros(ComplexF64, 9, 9, 1, grids.θ.N * grids.θ.f_quad, grids.ζ.N * grids.ζ.f_quad)
end


"""
    local_matrix_size(grids::FSSGridsT)

Computes the size of the local I, W matrices needed at each grid point.
"""
function local_matrix_size(grids::FSSGridsT)

    return zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.N * grids.θ.f_quad, grids.ζ.N * grids.ζ.f_quad)
end


"""
    local_matrix_size(grids::FFSGridsT)

Computes the size of the local I, W matrices needed at each grid point.
"""
function local_matrix_size(grids::FFSGridsT)

    return zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, grids.ζ.N * grids.ζ.f_quad)
end

"""
    local_matrix_size(grids::FFFGridsT)

Computes the size of the local I, W matrices needed at each grid point.
"""
function local_matrix_size(grids::FFFGridsT)

    return zeros(ComplexF64, 9, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)
end



"""
    init_bases_function(grids::FSSGridsT)

Initialises the empty structure that contains the trial and test functions.
4's are for number of Hermite elements when using FEM, 9 is for total first and second derivatives.
"""
function init_basis_function(grids::FSSGridsT)

    #order here is chosen so integration is done on contigeous block
    return zeros(ComplexF64, 4, 9, grids.r.gp)   
end



"""
    init_bases_function(grids::FFFGridsT)

Initialises the empty structure that contains the trial and test functions.
4's are for number of Hermite elements when using FEM, 9 is for total first and second derivatives.
"""
function init_basis_function(grids::FFFGridsT)

    #order here is chosen so integration is done on contigeous block
    return zeros(ComplexF64, 4, 4, 4, 9, grids.r.gp, grids.θ.gp, grids.ζ.gp)
end


"""
    init_bases_function(grids::FFSGridsT)

Initialises the empty structure that contains the trial and test functions.
4's are for number of Hermite elements when using FEM, 9 is for total first and second derivatives.
"""
function init_basis_function(grids::FFSGridsT)

    #order here is chosen so integration is done on contigeous block
    return zeros(ComplexF64, 4, 4, 9, grids.r.gp, grids.θ.gp) 
end

