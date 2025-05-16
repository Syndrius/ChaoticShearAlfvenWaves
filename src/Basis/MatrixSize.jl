"""
    matrix_size(grids::FSSGridsT)

Computes the dimension of the matrix.
"""
function matrix_size(grids::ContGridsT)

    return grids.x2.N * grids.x3.N
end


"""
    matrix_size(grids::FSSGridsT)

Computes the dimension of the matrix.
"""
function matrix_size(grids::FSSGridsT)

    return 2 * grids.x1.N * grids.x2.N * grids.x3.N
end


"""
    matrix_size(grids::FFSGridsT)

Computes the dimension of the matrix.
"""
function matrix_size(grids::FFSGridsT)

    return 4 * grids.x1.N * grids.x2.N * grids.x3.N
end


"""
    matrix_size(grids::FFFGridsT)

Computes the dimension of the matrix.
"""
function matrix_size(grids::FFFGridsT)

    return 8 * grids.x1.N * grids.x2.N * grids.x3.N
end

"""
    local_matrix_size(grids::FSSGridsT)

Computes the size of the local I, W matrices needed at each grid point.
"""
function local_matrix_size(grids::ContGridsT)

    return zeros(ComplexF64, 9, 9, 1, grids.x2.N * grids.x2.f_quad, grids.x3.N * grids.x3.f_quad)
end


"""
    local_matrix_size(grids::FSSGridsT)

Computes the size of the local I, W matrices needed at each grid point.
"""
function local_matrix_size(grids::FSSGridsT)

    return zeros(ComplexF64, 9, 9, grids.x1.gp, grids.x2.N * grids.x2.f_quad, grids.x3.N * grids.x3.f_quad)
end


"""
    local_matrix_size(grids::FFSGridsT)

Computes the size of the local I, W matrices needed at each grid point.
"""
function local_matrix_size(grids::FFSGridsT)

    return zeros(ComplexF64, 9, 9, grids.x1.gp, grids.x2.gp, grids.x3.N * grids.x3.f_quad)
end

"""
    local_matrix_size(grids::FFFGridsT)

Computes the size of the local I, W matrices needed at each grid point.
"""
function local_matrix_size(grids::FFFGridsT)

    return zeros(ComplexF64, 9, 9, grids.x1.gp, grids.x2.gp, grids.x3.gp)
end



"""
    init_bases_function(grids::FSSGridsT)

Initialises the empty structure that contains the trial and test functions.
4's are for number of Hermite elements when using FEM, 9 is for total first and second derivatives.
"""
function init_basis_function(grids::FSSGridsT)

    #order here is chosen so integration is done on contigeous block
    return zeros(ComplexF64, 4, 9, grids.x1.gp)   
end



"""
    init_bases_function(grids::FFFGridsT)

Initialises the empty structure that contains the trial and test functions.
4's are for number of Hermite elements when using FEM, 9 is for total first and second derivatives.
"""
function init_basis_function(grids::FFFGridsT)

    #order here is chosen so integration is done on contigeous block
    return zeros(ComplexF64, 4, 4, 4, 9, grids.x1.gp, grids.x2.gp, grids.x3.gp)
end


"""
    init_bases_function(grids::FFSGridsT)

Initialises the empty structure that contains the trial and test functions.
4's are for number of Hermite elements when using FEM, 9 is for total first and second derivatives.
"""
function init_basis_function(grids::FFSGridsT)

    #order here is chosen so integration is done on contigeous block
    return zeros(ComplexF64, 4, 4, 9, grids.x1.gp, grids.x2.gp) 
end

