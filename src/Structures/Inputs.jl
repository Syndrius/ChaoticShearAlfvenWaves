#we may want to predefine these fkn sizes given we use them all over the place.
#=
"""
Struct storing Tempory Matrices used for memory efficient weak form.
"""
struct TMT
    C :: mat39
    D :: mat33
    jpar :: vec1
    F :: vec9
    Γ :: mat33
    dΓ :: mat333
    K :: vec6
    function TMT()
        #these should probably be static arrays for fk sake.
        new(zeros(mat39),
            zeros(mat33),
            zeros(vec1),
            zeros(vec9),
            zeros(mat33),
            zeros(mat333),
            zeros(vec6)
           )
    end
end
=#

struct TMT
    C :: Array{Float64, 2}
    D :: Array{Float64, 2}
    jpar :: Array{Float64, 1}
    F :: Array{Float64}
    Γ :: Array{Float64, 2}
    dΓ :: Array{Float64, 3}
    K :: Array{Float64}
    function TMT()
        new(zeros(3, 9), zeros(3, 3), zeros(1), zeros(9), zeros(3, 3), zeros(3, 3, 3), zeros(6))
    end
end

#unsure what else this will need yet
#qfm version of this will need more stuff naturally!
struct WeakFormInputsT
    B :: BFieldT
    met :: MetT
    tm :: TMT
    function WeakFormInputsT()
        new(BFieldT(), MetT(), TMT())
    end
end
