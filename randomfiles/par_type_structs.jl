
#testing if we can change our key strutcs to be parametric types of other struct.

abstract type BFieldT end
abstract type GeoT end
abstract type GridsT end
#%%

struct Prob{BT <: BFieldT, GeT <: GeoT, GrT <: GridsT}
    B :: BT
    Geo :: GeT
    Grids :: GrT
end
#%%

struct BaseBT <: BFieldT
    q :: Function
    a ::Float64
end

struct BaseGridsT <: GridsT
    N :: Int64
end

struct BaseGeoT <: GeoT
    met :: Int64
end
#%%


con_geo = BaseGeoT(4)
con_grids = BaseGridsT(1)
con_B = BaseBT(sin, 3.8)

con_prob = Prob(con_B, con_geo, con_grids)

#%%
#ok yeah this could be real nice actually.
function md_test(prob::Prob{BaseBT, <:GeoT, <:GridsT})

    display(prob.B.q(0.2))
end

typeof(con_prob)
md_test(con_prob)
