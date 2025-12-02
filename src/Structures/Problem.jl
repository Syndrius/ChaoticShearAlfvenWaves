"""

Main struct used for computing the weak form.
This struct is initialised in WeakForm.

### Fields
- fields::FieldsT - Stores the fields used.
- geo::GeometryT - Stores the geometry of the problem.
- flr::FLRT - Stores any finite Larmor radius effects.
"""
struct ProblemT 
    fields :: FieldsT
    geo :: GeometryT
    flr :: FLRT
end
