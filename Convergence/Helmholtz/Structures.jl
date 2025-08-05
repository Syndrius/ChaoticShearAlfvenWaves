
#basic structrue for the test cases

@kwdef struct TestProblemT <: ProblemT
    #doesn't actually need any
    #but this is used for Hermitian or not.
    flr :: FLRT = Structures.no_flr
    geo :: GeoParamsT
end


