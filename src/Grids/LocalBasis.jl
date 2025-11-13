#these are working now, save a little bit of allocations.
#could probably be combined into one thing!
#simple in place form , which is just passed on the apropriate shit being passed in.
#try generalise later!
function local_to_global!(x1::Array{Float64}, dx1::Array{Float64}, ind::Int64, ξ::Array{Float64}, grid::Array{Float64})

    dx1[1] = grid[ind+1] - grid[ind]

    #local coords are mapped from (-1, 1) to (0, 1)
    #then scaled by dx
    #finally, shiften to the proper place in the grid.
    @. x1 = dx1[1] * ((ξ+1)/2) + grid[ind]

    return dx1[1] / 2
end


function local_to_global!(x1::Array{Float64}, x2::Array{Float64}, dx::Array{Float64}, x1ind::Int64, x2ind::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen)

    dx[1] = x1grid[x1ind+1] - x1grid[x1ind]

    if x2ind == length(x2grid)
        dx[2] = 2π + x2grid[1] - x2grid[end]
    else
        dx[2] = x2grid[x2ind+1] - x2grid[x2ind]
    end

    #local coords are mapped from (-1, 1) to (0, 1)
    #then scaled by dx
    #finally, shiften to the proper place in the grid.
    @. x1 = dx[1] * ((ξx1+1)/2) + x1grid[x1ind]
    @. x2 = dx[2] * ((ξx2+1)/2) + x2grid[x2ind]

    return dx[1] * dx[2] / 4
end

#
function local_to_global!(x1::Array{Float64}, x2::Array{Float64}, x3::Array{Float64}, dx::Array{Float64}, x1ind::Int64, x2ind::Int64, x3ind::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, ξx3::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen, x3grid::StepRangeLen)

    dx[1] = x1grid[x1ind+1] - x1grid[x1ind]

    if x2ind == length(x2grid)
        dx[2] = 2π + x2grid[1] - x2grid[end]
    else
        dx[2] = x2grid[x2ind+1] - x2grid[x2ind]
    end

    if x3ind == length(x3grid)
        dx[3] = 2π + x3grid[1] - x3grid[end]
    else
        dx[3] = x3grid[x3ind+1] - x3grid[x3ind]
    end

    #local coords are mapped from (-1, 1) to (0, 1)
    #then scaled by dx
    #finally, shiften to the proper place in the grid.
    @. x1 = dx[1] * ((ξx1+1)/2) + x1grid[x1ind]
    @. x2 = dx[2] * ((ξx2+1)/2) + x2grid[x2ind]
    @. x3 = dx[3] * ((ξx3+1)/2) + x3grid[x3ind]

    return dx[1] * dx[2] * dx[3] / 8
end

############################################

#we will just need to figure out how to specify a arbitrary sized tuple, think Ntuple is the go.
#tuple here is weird af.
#types for 1d are a bit weird because the tuples don't really work.
#either have to force them to be tuples, or just go with different defs.
function local_to_global!(r::Tuple, dr::Array{Float64}, ind::CartesianIndex, ξ::Vector{Float64}, grids::AbstractArray{Float64})

    #should be able to work for all cases now.
    #here I think we assume periodicicty
    #then the constraint is based on which grid points are iterated through
    #this case does not generalise properly!
    for i in 1:1 #perhaps length of of dr?
        n = ind[i]
        #if n == length(grids[i])
        if n == length(grids)
            #ideally, something like this, but these are the instantiated grids.
            #more and more like these should just be stored as part of it!
            #dr[i] = grids.stop + grids[i][1] - grids[i][n]
            #dr[i] = 2π + grids[i][1] - grids[i][n]
            dr[i] = 2π + grids[1] - grids[n]
        else
            #dr[i] = grids[i][n+1] - grids[i][n]
            #need to properly put this grids in the same form for this to work!
            dr[i] = grids[n+1] - grids[n]
        end
        #jn = ind[i]
        #dr[i] = grids[i][n+1] - grids[i][n]
        #diff as grids is not a tuple for this case.
        #dr[i] = grids[n+1] - grids[n]

        
        #display(i)
        #display(r)
        @. r[i] = ((ξ[i]+1) / 2) * dr[i] + grids[n]
    end

    return dr[1] / 2

end

#more generic version, that should cater towards sm
function local_to_global!(r::NTuple{2, Array{Float64}}, dr::Array{Float64}, ind::CartesianIndex, ξ::NTuple{2, Array{Float64}}, grids::NTuple{2, Array{Float64}})
    display("this one!")

    #ξ will only be the length of fem cases
    #this will only work if the last dim in sm. i.e. we wont' be able to do fem sm fem.
    for i in 1:length(dr)
        n = ind[i]
        if n == length(grids[i])
            #dr[i] = grids.stop + grids[i][1] - grids[i][n]
            dr[i] = 2π + grids[i][1] - grids[i][n]
        else
            dr[i] = grids[i][n+1] - grids[i][n]
        end

        display(length.(r))
        display(length.(ξ))
        display(length.(grids))
        @. r[i] = ((ξ[i]+1) / 2) * dr[i] + grids[i][n]
    end
    jac = 1.0
    for i in 1:length(dr)
        jac *= dr[i]
    end
    return jac / (2^length(dr))
end

#these tuple types are cppled af.
#function local_to_global!(r::Tuple{Array{Float64}, Array{Float64}}, dr::Array{Float64}, ind::CartesianIndex, ξ::Tuple{Array{Float64}, Array{Float64}}, grids::Tuple{AbstractArray{Float64}, AbstractArray{Float64}})
#unsure if allowing for floats and ints is a good or bad idea...
function local_to_global!(r::Tuple{AbstractArray, AbstractArray}, dr::Array{Float64}, ind::CartesianIndex, ξ::Tuple{AbstractArray, AbstractArray}, grids::Tuple{AbstractArray{Float64}, AbstractArray{Float64}})

    #should be able to work for all cases now.
    #here I think we assume periodicicty
    #then the constraint is based on which grid points are iterated through
    #so this should be indep of SM basis...
    #i.e. nothing happens here for SM.
    #so now we are using dr to determine if we need to do the transformation
    #bit silly, but might work.
    for i in 1:length(dr) #perhaps length of of dr?
        n = ind[i]
        if n == length(grids[i])
            #dr[i] = grids.stop + grids[i][1] - grids[i][n]
            dr[i] = 2π + grids[i][1] - grids[i][n]
        else
            dr[i] = grids[i][n+1] - grids[i][n]
        end

        @. r[i] = ((ξ[i]+1) / 2) * dr[i] + grids[i][n]
    end
    jac = 1.0
    for i in 1:length(dr)
        jac *= dr[i]
    end
    return jac / (2^length(dr))
end

function local_to_global!(r::Tuple{Array{Float64}, Array{Float64}, Array{Float64}}, dr::Array{Float64}, ind::CartesianIndex, ξ::Tuple{Array{Float64}, Array{Float64}, Array{Float64}}, grids::Tuple{AbstractArray{Float64}, AbstractArray{Float64}, AbstractArray{Float64}}) #where (A::Array{Float64})

    #should be able to work for all cases now.
    #here I think we assume periodicicty
    #then the constraint is based on which grid points are iterated through
    for i in 1:3 #perhaps length of of dr?
        n = ind[i]
        if n == length(grids[i])
            #dr[i] = grids.stop + grids[i][1] - grids[i][n]
            dr[i] = 2π + grids[i][1] - grids[i][n]
        else
            dr[i] = grids[i][n+1] - grids[i][n]
        end
        #n = ind[i]
        #dr[i] = grids[i][n+1] - grids[i][n]

        #display(r)
        @. r[i] = ((ξ[i]+1) / 2) * dr[i] + grids[i][n]
    end

    return dr[1] * dr[2] * dr[3] / 8

end

function local_to_global!(r::Tuple{Array{Float64}, Array{Float64}}, dr::Tuple{Float64, Float64}, ind::CartesianIndex, ξ1::Array{Float64}, ξ2::Array{Float64}, x1grid::AbstractArray{Float64}, x2grid)
    n1 = ind[1]
    n2 = ind[2]
    #this will need the periodicity constraint unfor.
    dr[1] = x1grid[n1+1] - x1grid[n1]
    dr[2] = x2grid[n2+1] - x2grid[n2]

    @. r[1] = ((ξ1+1) / 2) * dr[1] + x1grid[n1]
    @. r[2] = ((ξ2+1) / 2) * dr[2] + x2grid[n2]
end


#this is perhaps part of basis?
#either way this needs a huge amount of work.
#this is cooked af.
function local_to_global(ind::CartesianIndex, ξ::Array{Float64}, grid::AbstractArray{Float64})
    #think the benefit of Cartesian index here is we should be able to do each dim with the same function.
    node = ind[1]
    dx1 = grid[node+1] - grid[node]

    #first we map the (-1, 1) global coords to (0, 1)
    mp = @. (ξ+1) / 2

    #then we sale by the width of grid points
    mp = mp .* dx1

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    x1global = mp .+ grid[node]

    
    return x1global, dx1
end 

"""
    local_to_global(node::Int64, ξ::Array{Float64}, grid::Array{Float64})

Converts the local grid where the finite element basis is defined to the global coordinates.
"""
function local_to_global(node::Int64, ξ::Array{Float64}, grid::Array{Float64})

    #TODO
    #this could be made non-allocating

    #should this be in basis??
    dx1 = grid[node+1] - grid[node]

    #first we map the (-1, 1) global coords to (0, 1)
    mp = @. (ξ+1) / 2

    #then we sale by the width of grid points
    mp = mp .* dx1

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    x1global = mp .+ grid[node]

    
    return x1global, dx1
end 




"""
    local_to_global(x1node::Int64, x2node::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen)

Converts the local grid where the finite element basis is defined to the global coordinates.
"""
function local_to_global(x1node::Int64, x2node::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen)

    #handles periodicity
    if x2node == length(x2grid)
        dx2 = 2π + x2grid[1] - x2grid[x2node]
    else
        dx2 = x2grid[x2node+1] - x2grid[x2node]
    end


    dx1 = x1grid[x1node+1] - x1grid[x1node]
    

    #first we map the (-1, 1) global coords to (0, 1)
    mpx1 = @. (ξx1+1) / 2
    mpx2 = @. (ξx2+1) / 2

    #then we sale by the width of grid points
    mpx1 = mpx1 .* dx1
    mpx2 = mpx2 .* dx2

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    x1global = mpx1 .+ x1grid[x1node]
    x2global = mpx2 .+ x2grid[x2node]

    
    return x1global, x2global, dx1, dx2
end 


"""
    local_to_global(x1node::Int64, x2node::Int64, x3node::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, ξx3::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen, x3grid::StepRangeLen)

Converts the local grid where the finite element basis is defined to the global coordinates.
"""
function local_to_global(x1node::Int64, x2node::Int64, x3node::Int64, ξx1::Array{Float64}, ξx2::Array{Float64}, ξx3::Array{Float64}, x1grid::Array{Float64}, x2grid::StepRangeLen, x3grid::StepRangeLen)

    #handles periodicity
    if x2node == length(x2grid)
        dx2 = 2π + x2grid[1] - x2grid[x2node]
    else
        dx2 = x2grid[x2node+1] - x2grid[x2node]
    end

    if x3node == length(x3grid)
        dx3 = 2π + x3grid[1] - x3grid[x3node]
    else
        dx3 = x3grid[x3node+1] - x3grid[x3node]
    end

    dx1 = x1grid[x1node+1] - x1grid[x1node]
    
    #first we map the (-1, 1) global coords to (0, 1)
    mpx1 = @. (ξx1+1) / 2
    mpx2 = @. (ξx2+1) / 2
    mpx3 = @. (ξx3+1) / 2

    #then we sale by the width of grid points
    mpx1 = mpx1 .* dx1
    mpx2 = mpx2 .* dx2
    mpx3 = mpx3 .* dx3

    #finally we add this to our grid, so this gives us go points between our global grid point i and i+1
    x1global = mpx1 .+ x1grid[x1node]
    x2global = mpx2 .+ x2grid[x2node]
    x3global = mpx3 .+ x3grid[x3node]

    
    return x1global, x2global, x3global, dx1, dx2, dx3
end 
