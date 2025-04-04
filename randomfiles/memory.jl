

#so precompiling is probbaly going to make a fkn huge difference
#compilation is using 0.25 of a Gb (ish).
#if this is true for each proc, that is a huge amount of memory!

#precom for MIDParallel is using ~0.6 Gb.
#so compilation could be using almost a Gb per proc...
#that will make a huge difference to our memory usage...
#but then again this memory should be reused so who cares???

#so to summarise, using structs causes allocations,
#using views/submatrices causes allocations.
#not sure we can fix any of this.
#guess we should consider solve independently, and consider 
#precompilation.
#still need to investigate further what is using all the damn memory...
#may need to investigate the parallel case in serial, to see what is going on with the petsc calls.


@time using MIDParallel
@time using MID
using BenchmarkTools

@time function run_MID()
    


    Nr = 40;
    Nθ = 6

    geo = GeoParamsT(R0=10.0)

    prob = init_problem(q=Axel_q, geo=geo); 

    rgrid = rfem_grid(N=Nr)
    θgrid = afem_grid(N=Nθ, pf=2)
    ζgrid = asm_grid(start=-2, N=1)
    #grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
    grids = init_grids(rgrid, θgrid, ζgrid)

    evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true);

end

@code_warntype run_MID();



Nr = 40;
Nθ = 6

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo); 

rgrid = rfem_grid(N=Nr)
θgrid = afem_grid(N=Nθ, pf=2)
ζgrid = asm_grid(start=-2, N=1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);
grids = init_grids(rgrid, θgrid, ζgrid)

@allocated W, Imat = construct(prob, grids)


#allocated for full process 1200981672 ~1.2Gb
#for construct 455466848 ~0.45Gb
#changing the met and B contruct is now 319784512 ~0.31Gb.
#still room for improvement here, i.e. there is still a fair bit of allocation going on.
#for solve 88592688 ~0.08Gb
#post_process 420328824 ~0.4Gb


@allocated evals, efuncs = full_spectrum_solve(Wmat=W, Imat=Imat, ideal=true)

@allocated evals, ϕ, ϕft = MID.Spectrum.post_process(evals, efuncs, grids, prob.geo, false)

@allocated MID.Spectrum.matrix_size(grids)

@allocated ϕ_g, ϕft_g = MID.PostProcessing.allocate_phi_arrays(grids, length(evals.ω), deriv=false)


@time construct(prob, grids);

@time rgrid, θgrid, ζgrid = inst_grids(grids)

@time MID.Structures.inst_grid(grids.θ)


@benchmark MID.Structures.inst_grid(grids.θ)

@time B = MID.MagneticField.BFieldT()
@time met = MID.Geometry.MetT()

9*8 + 9*8 + 27*8 + 27*8 + 8 * 3*8

@time ξr, wgr = MID.Spectrum.gausslegendre(grids.r.gp) 
ξθ, wgθ = MID.Spectrum.gausslegendre(grids.θ.gp)

@time S = MID.Basis.hermite_basis(ξr, ξθ);

@time Φ = MID.Basis.init_basis_function(grids);

using FFTW

@time Wmat = MID.Basis.local_matrix_size(grids);

@time p = plan_fft!(Wmat, [5])


@time r, θ, dr, dθ = MID.Basis.local_to_global(1, 1, ξr, ξθ, rgrid, θgrid)

@time tm=MID.WeakForm.TM();

@time Imat = MID.Basis.local_matrix_size(grids);
@time Wmat = MID.Basis.local_matrix_size(grids);

@benchmark MID.WeakForm.W_and_I!(Wmat, Imat, met, B, prob, r, θ, ζgrid, tm)

@benchmark prob.compute_met(met, r[1], θ[1], ζgrid[1], prob.geo.R0)

@benchmark Δp = r[1]/(4*prob.geo.R0)

@benchmark MID.MagneticField.compute_B!(B, met, prob.q, prob.isl, r[1], θ[2], ζgrid[3])

@benchmark MID.WeakForm.compute_D!(met, B, tm.D)

@benchmark @views MID.WeakForm.compute_W!(Wmat[:, :, 1, 2, 3], met, B, 1.0, 0.0, tm)

@code_warntype @views MID.WeakForm.compute_W!(Wmat[:, :, 1, 2, 3], met, B, 1.0, 0.0, tm)


@benchmark @views MID.WeakForm.compute_I!(Imat[:, :, 1, 2, 3], met, B, 1.0, prob.flr, tm.D, tm.F)


@benchmark @views MID.WeakForm.compute_F!(met, B, tm.F)

@benchmark @views MID.WeakForm.Tl!(Wmat[:, :, 1, 2, 3], met, B, tm.C, tm.D, tm.T)

struct testStruct
    a :: Float64
end

@kwdef struct kStruct
    b :: Float64
end

ms = testStruct(0.65)
mks = kStruct(b=0.65)
function test1(s::testStruct)

    

    return s.a + 6
end

function test2(s::kStruct)

    return s.b+6
end

function test3(c::Float64)

    return c + 6
end


@benchmark test1(ms)
@benchmark test2(mks)
@benchmark test3(mks.b)
@benchmark test3(ms.a)



function test_toroidal_metric!(met::MID.Geometry.MetT, r::Float64, θ::Float64, ζ::Float64, prob::ProblemT)#R0::Float64)#

    """
    
    Δp = r/(4*R0)

    ϵ = r/R0


    η = 1/2 * (ϵ+Δp)

    met.gl[1, 1] = 0.1 * Δp + 4 * r
    met.gl[1, 2] = 0.1 * cos(θ)
    met.J = η * 2
    met.dJ[1] = 2*η
end
    
    """
    #R0 = prob.geo.R0

    Δp = r/(4*prob.geo.R0)
    Δpp = 1/(4*prob.geo.R0)

    ϵ = r/prob.geo.R0

    η = 1/2 * (ϵ+Δp)

    ηp = 1/2* (1/prob.geo.R0 + Δpp)

    #met.J = r * prob.geo.R0 * (1+2*ϵ*cos(θ))

    met.gl[1, 1] = 1-2*Δp * cos(θ)
    met.gl[1, 2] = r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gl[2, 1] = r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gl[2, 2] = r^2*(1+4*η*cos(θ) + 4*η^2)
    met.gl[3, 3] = prob.geo.R0^2*(1+2*ϵ*cos(θ))


    met.gu[1, 1] = 1+2*Δp * cos(θ)
    met.gu[1, 2] = -1/r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gu[2, 1] = -1/r*(ϵ + Δp + r*Δpp) * sin(θ)
    met.gu[2, 2] = 1/r^2*(1-2*(ϵ+Δp)*cos(θ))
    met.gu[3, 3] = 1/prob.geo.R0^2*(1-2*ϵ*cos(θ))

    
    met.dJ[1] = prob.geo.R0 + 4*r * cos(θ)
    met.dJ[2] = -2 * r * prob.geo.R0*ϵ * sin(θ)

    #first two indicies give metric element, while third is derivative,
    #eg [1, 2, 3] is ∂g_{12}/∂ζ
    met.dgl[1, 1, 1] = -2*Δpp * cos(θ)
    met.dgl[1, 1, 2] = 2*Δp * sin(θ)

    met.dgl[1, 2, 1] = ((ϵ + Δp + r*Δpp) + r*(1/prob.geo.R0 + 2*Δpp)) * sin(θ)
    met.dgl[1, 2, 2] = r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgl[2, 1, 1] = ((ϵ + Δp + r*Δpp) + r*(1/prob.geo.R0 + 2*Δpp)) * sin(θ)
    met.dgl[2, 1, 2] = r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgl[2, 2, 1] = 2*r*(1+4*η*cos(θ) + 4*η^2) + r^2 * (4*ηp*cos(θ) + 8*η * ηp)
    met.dgl[2, 2, 2] = -r^2*(4*η*sin(θ))

    met.dgl[3, 3, 1] = 2*prob.geo.R0*cos(θ)
    met.dgl[3, 3, 2] = -2*prob.geo.R0^2*ϵ*sin(θ)

    met.dgu[1, 1, 1] = 2*Δpp * cos(θ)
    met.dgu[1, 1, 2] = -2*Δp * sin(θ)

    met.dgu[1, 2, 1] = (1/r^2 * (ϵ + Δp + r*Δpp) - 1/r*(1/prob.geo.R0 + 2*Δpp)) * sin(θ)
    met.dgu[1, 2, 2] = -1/r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgu[2, 1, 1] = (1/r^2 * (ϵ + Δp + r*Δpp) - 1/r*(1/prob.geo.R0 + 2*Δpp)) * sin(θ)
    met.dgu[2, 1, 2] = -1/r*(ϵ + Δp + r*Δpp) * cos(θ)

    met.dgu[2, 2, 1] = -2/r^3 * (1-2*(ϵ+Δp)*cos(θ)) + 1/r^2 * (-2*(1/prob.geo.R0+Δpp)*cos(θ))
    met.dgu[2, 2, 2] = 2/r^2*(ϵ+Δp)*sin(θ)

    met.dgu[3, 3, 1] = -2*cos(θ)/prob.geo.R0^3
    met.dgu[3, 3, 2] = 2*ϵ*sin(θ)/prob.geo.R0^2


end
#"""

#dJ is not type stable for some reason????
#perhaps becuase we havent specified the 1? -> fixed that issue


#this may not be worth the effort tbh.
typeof(@view r[1])
typeof(@views(r[1]))

@benchmark test_toroidal_metric!(met,  getindex(r, 1), getindex(θ, 1), getindex(ζgrid, 1), 10.0)#prob.geo.R0)

R0 = 10.0

tmp = geo.R0

#@benchmark test_toroidal_metric!(met, r[1], θ[1], ζgrid[1], geo.R0)
@code_warntype test_toroidal_metric!(met, r[1], θ[1], ζgrid[1], prob)
@benchmark test_toroidal_metric!(met, 0.1, 0.2, 0.3, prob)

@profview test_toroidal_metric!(met, 0.1, 0.2, 0.3, prob)


@benchmark test_toroidal_metric!(met, 0.2, 0.3, 0.4, 10.0)
@code_warntype test_toroidal_metric!(met, 0.2, 0.3, 0.4, @views prob.geo.R0)




@benchmark MID.Basis.create_local_basis!(Φ, S, grids.θ.pf, 2, dr, dθ)
@benchmark MID.Basis.create_local_basis!(Φ, S, 3, 2, dr, dθ)
@code_warntype MID.Basis.create_local_basis!(Φ, S, grids.θ.pf, 2, dr, dθ)


@benchmark right_ind = MID.Indexing.grid_to_index(1, 5, 3, 1, 3, grids)

testr = 2
testθ = 1
trialr = 3
trialθ = 1
nind=1
jac = 0.1/2

#@benchmark Wsum = @views MID.Integration.gauss_integrate(Φ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], Wmat[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)

@benchmark Wsum = @views MID.Integration.gauss_integrate(Φ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], Wmat[:, :, :, :, nind], wgr, wgθ, jac, 4, 4)

display(grids.r.gp)


function gauss_integrate_test(Ψ::Array{ComplexF64, 3}, Φ::Array{ComplexF64, 3}, mat::Array{ComplexF64, 4}, wgr::Array{Float64}, wgθ::Array{Float64}, jac::Float64, rgp::Int64, θgp::Int64)

    res = 0.0 + 0.0im
    for l in 1:θgp, k in 1:rgp

        scale = wgr[k] * wgθ[l] * jac

        for j in 1:9, i in 1:9
            res += @inbounds Ψ[i, k, l] * mat[i, j, k, l] * Φ[j, k, l] * scale
        end
    end

    return res
end

@benchmark @views Wsum = gauss_integrate_test(Φ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], Wmat[:, :, :, :, nind], wgr, wgθ, jac, 4, 4)


@code_warntype gauss_integrate_test(view(Φ, testr, testθ, :, :, :), view(Φ, trialr, trialθ, :, :, :), view(Wmat, :, :, :, :, nind))#, wgr, wgθ, jac, 4, 

@benchmark gauss_integrate_test($view(Φ, testr, testθ, :, :, :), $view(Φ, trialr, trialθ, :, :, :), $view(Wmat, :, :, :, :, nind))


#why the fk is this allocating anything????
@views function gauss_integrate_test(Ψ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, Φ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, mat::SubArray{ComplexF64, 4, Array{ComplexF64, 5}})#, wgr::Array{Float64}, wgθ::Array{Float64}, jac::Float64, rgp::Int64, θgp::Int64)

    #res = 0.0 + 0.0im
    #k = 1
    #l = 1

    return 2.0 #Φ[2, 1, 3]# + mat[2, 1, 3, 1]

       

    #return res
end


function integrate(a::Array{ComplexF64, 4})

    return a[1, 2, 1, 2]
end

@benchmark @views integrate(Wmat[1, :, :, :, :])

@benchmark integrate(A)

A = zeros(ComplexF64, 3, 4, 2, 3)
B = zeros(3, 4, 2)

#not useful I think...
using Tullio


#allocated twice as much and takes ~3 times longer...
#not ideal.
function integrate_tul(Ψ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, Φ::SubArray{ComplexF64, 3, Array{ComplexF64, 5}}, mat::SubArray{ComplexF64, 4, Array{ComplexF64, 5}}, wgr::Array{Float64}, wgθ::Array{Float64}, jac::Float64, rgp::Int64, θgp::Int64)

    @tensoropt res = Ψ[i, k, l] * mat[i, j, k, l] * Φ[j, k, l] * wgr[k] * wgθ[l] * jac

    return res



end


@benchmark @views Wsum = integrate_tul(Φ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], Wmat[:, :, :, :, nind], wgr, wgθ, jac, 4, 4)
