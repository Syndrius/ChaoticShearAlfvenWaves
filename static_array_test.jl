
#just trying to get a struct full of static arrays to work!

using StaticArrays

#%%
struct test2{L, T}
    x :: SArray{L, T}
end

#%%
mytest = test{10}


data = LinRange(0, 1, 10)


mytest = test{2}#, data}

data = SArray{2}([1 2])

data = @SArray [1, 2]

mytest.x = data

mytest = test{2}(@SArray [1.0, 2.0])

mytest = test2{2, Float64}
#so this doesn't work even if we could get the damn structs to work...
isbitstype(mytest)
