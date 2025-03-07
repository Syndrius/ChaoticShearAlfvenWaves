
#no idea what to call this damn thing!
#not even sure we want this tbh!
struct CoefficientsT
    nv :: Array{Float64, 1} #this is actually just a float. but we want to mutate it.
    rcos :: Array{Float64, 1}
    θcos :: Array{Float64, 1}
    rsin :: Array{Float64, 1}
    θsin :: Array{Float64, 1}
    function CoefficientsT(N)
        new(zeros(1), zeros(N), zeros(N), zeros(N), zeros(N))
    end
end


#redoing these functions but less terribly.

#unsure if action is the best name for this but whatever.
#this is now much much better.
#gets better the larger the functions are 
function action2(p::Int64, q::Int64, prob::ProblemT, met::MetT, B::BFieldT)#, sguess::Float64, MM::Int64, M::Int64, N::Int64)

    MM = 4 #this has to be 2 * nfft_multiplier, which is defaulted to 2 in rfft functions.
    #but this at 4 and pqNtor the function is v quick.
    #Ntor = 10
    N = 8
    M = 24 #only used in the very last step, so something to do with α?
    sguess = 1.1
    #these are some of the most unclear peices of garbage going around.
    qN = q * N
    fM = MM * N
    qfM = q * fM
    Nfft= qfM

    Ntor = qN + 1 #actual size of our fft arrays.

    #no idea what these are either!
    dζ = 2π / fM
    dθ = 2π / qfM

    

    nlist = range(0, qN)
    #don't think zeta is the toroidal component.
    ζ = range(0, Nfft-1) .* dζ

    nzq = nlist .* ζ' ./ q

    #action array.
    #minus two for the two sin constant terms
    #which we don't solve for.
    #our definintion doesn't matter!
    #it just creates one when it wants!
    #δS = zeros(4*Ntor+1 - 2)

    #array that stores the all the coefficients for root finding
    #-2 because we don't include the sin 0 coeffs.
    x0 = zeros(4*Ntor+1-2)

    #x0 = 

    #we take our r, θ original toroidal coordinates as an arbitrary fourier expansion in terms of zeta to denote our surface
    #r(ζ) = r_0^c ∑_{n=1}^qN (r_n^c cos(nζ/q) + r_n^s sin(nζ/q))
    #θ(ζ) = θ_0^c ∑_{n=1}^qN (θ_n^c cos(nζ/q) + θ_n^s sin(nζ/q))

    #in a very ideal world, a lot of this would be defined outside of this function, especially if we are constructing many coords.

    #these store the minimised values for the coefficeints at each iteration
    #need to work our what the damn iteration is.
    rcosarr = zeros(fM , Ntor)
    rsinarr = zeros(fM , Ntor)
    θcosarr = zeros(fM , Ntor)
    θsinarr = zeros(fM , Ntor)
    nvarr = zeros(fM)

    

    nfft_mulitplier = 2

    #unsure about this tbh.
    size_of_fft_array = nfft_mulitplier * (Ntor - 1)

    #need to add 1 here so irft doesn't screem.
    ift_r1D = zeros(ComplexF64, size_of_fft_array+1)
    ift_θ1D = zeros(ComplexF64, size_of_fft_array+1)

    #predefine fourier transform arrays and plan for efficiency
    #not actually sure on the shape
    #so input rfft needs to be size_of_fft +1 and
    #output needs to be 2*size_of_fft_array.
    #not really clear what is going on here with the real fft.
    r1D = zeros(2*size_of_fft_array)
    θ1D = zeros(2*size_of_fft_array)

    #arrays for storing the magnetic field after it is computed one coord at a time.
    Br = zeros(2*size_of_fft_array)
    Bθ = zeros(2*size_of_fft_array)
    Bζ = zeros(2*size_of_fft_array)

    #size of this is based on a single example. not true in general.
    ft_r1D = zeros(ComplexF64, size_of_fft_array + 1)
    ft_θ1D = zeros(ComplexF64, size_of_fft_array + 1)


    rdot_fft_cos = zeros(size_of_fft_array+1)
    rdot_fft_sin = zeros(size_of_fft_array+1)

    θdot_fft_cos = zeros(size_of_fft_array+1)
    θdot_fft_sin = zeros(size_of_fft_array+1)

    #this can't be done in place apparently!
    #maybe something to do with the padding issue we will face.
    irft_1D_p = plan_irfft(ift_r1D, 2*size_of_fft_array)
    #display(length(ft_r1D))

    #display(length(Br))
    rft_1D_p = plan_rfft(Br)

    #srs what is this loop
    #each of these gives a solution, but this doesn't seem to be hte number of surfaces?
    #I guess this might be number of surfaces times number of fourier harmonics in alpha?
    #no, this is for a single surface.
    #no fkn idea what this is looping over.
    #struct for storing the coefficients
    #CT sttruct seems to cause the solver to freak out,
    #spends ages checking for singularity.
    CT = CoefficientsT(qN+1)
    for jpq in 0:fM-1
        #maybe this is supposed to be alpha?
        #or area? fk knows.
        a = jpq * dθ
        if jpq == 0
            
            #initial guesses for the alg.
            CT.rcos[1] = sguess
            #everything else is zero.
            rcos0 = zeros(qN+1)
            rsin0 = zeros(qN+1)
            θcos0 = zeros(qN+1)
            θsin0 = zeros(qN+1)

            rcos0[1] = sguess
            θcos0[1] = 0

            nv0 = 0


        else
            #otherwise use the last solution as a starting point.
            CT.rcos .= rcosarr[jpq, :]
            CT.rsin .= rsinarr[jpq, :]
            CT.θcos .= θcosarr[jpq, :]
            CT.θsin .= θsinarr[jpq, :]
            CT.nv[1] = nvarr[jpq]
            rcos0 = rcosarr[jpq, :]
            rsin0 = rsinarr[jpq, :]
            θcos0 = θcosarr[jpq, :]
            θsin0 = θsinarr[jpq, :]
            nv0 = nvarr[jpq]

            θcos0[1]  += dθ 
            CT.θcos[1] += dθ
        end
        #combine the values into a single array for the root findning alg's argument
        #x0[1:Ntor] = rcos0
        #x0[Ntor+1:2*Ntor-1] = rsin0[2:end]
        #x0[2*Ntor+1-1:3*Ntor-1] = θcos0
        #x0[3*Ntor+1-1:4*Ntor-2] = θsin0[2:end]
        #x0[end] = nv0
        
        #x0 = pack_dof(nv0, rcos0, θsin0, rsin0, θcos0)
        pack_coeffs!(x0, CT, Ntor)

        #x0 .= 0.0
        #x0[3] = 8.0
        #first element cannot be zero!
        #x0[1] = 2.0
        #x0[2] = 3.0
        #x0[3] = 4.0
        #x0[end-1] = 7.0
        #x0[3*Ntor] = 7.0


        #display(ff)
        #display(x0)
        #display(p)
        #display(q)
        #display(a)
        #display(qN)
        #display(ζ)
        #display(nlist)
        #display(prob)

        #reduce args of function for root finding
        #this needs fkn heaps of args.
        #may want to create a struct for this. this is getting ridic.
        
        #so the in place version doesn't seem to work???
        #in place is not the problemo!
        #we had forgotten to return the arg lol.
        #think it is just not working at all lol
        ag!(δS, x) = action_grad!(δS, x, CT, p, q, a, Ntor, r1D, θ1D, ζ, irft_1D_p, prob, met, B, Br, Bθ, Bζ, ift_r1D, ift_θ1D, nlist, rdot_fft_cos, rdot_fft_sin, θdot_fft_cos, θdot_fft_sin, ft_r1D, ft_θ1D, rft_1D_p)
        #ag!(x) = action_grad!(x, p, q, a, Ntor, r1D, θ1D, ζ, irft_1D_p, prob, met, B, Br, Bθ, Bζ, ft_r1D, ft_θ1D, nlist)
        #ag!(δS, x) = action_gradient(δS, x, p, q, a, qN, ζ, nlist, prob)

        #so δS is different somewhere
        #need to figure out what is going on
        #think we have created the simplest problem now.
        #ag!(δS, x0)
        #δS = ag!(x0)
        #display(δS)#[1:5])
        
        sol = nlsolve(ag!, x0)

        #display(sol.zero)
        #nv, rcos, tsin, rsin, tcos = unpack_dof(sol.zero, qN)

        unpack_coeffs!(sol.zero, CT, Ntor)
        #rcos = sol.zero[1:Ntor]
        #rsin = [[0] ; sol.zero[Ntor+1:2*Ntor-1]]
        #tcos = sol.zero[2*Ntor+1-1:3*Ntor-1]
        #tsin = [[0] ; sol.zero[3*Ntor+1-1:4*Ntor - 2]]

        #rcosarr[jpq+1, :] = sol.zero[1:Ntor]
        #rsinarr[jpq+1, :] = sol.zero[Ntor+1:2*Ntor]
        #θcosarr[jpq+1, :] = sol.zero[2*Ntor+1:3*Ntor]
        #θsinarr[jpq+1, :] = sol.zero[3*Ntor+1:4*Ntor]
        rcosarr[jpq+1, :] = CT.rcos
        rsinarr[jpq+1, :] = CT.rsin
        θcosarr[jpq+1, :] = CT.θcos
        θsinarr[jpq+1, :] = CT.θsin
        nvarr[jpq+1] = CT.nv[1]
        
    end

    # I think this should be in a different function
    #this is not the same as solving the action
    #this is converting the minimised surface to the straight field line coords I think.
    r = irfft1D(rcosarr, rsinarr)
    ζ = LinRange(0, 2*q*π, size(r)[end]+1)[1:end-1]

    
    θcosarr[:, 1] .= 0.0

    

    θ = irfft1D(θcosarr, θsinarr)

    #display(t)

    r2D_alpha = zeros((qfM, Nfft))
    θ2D_alpha = zeros((qfM, Nfft))

    #this indexing is cooked af.
    for i in 0:q-1
        idx = mod(p*i, q)

        #display(i)

        #display(r[:, 1: i * pqNtor * MM])
        #display(r[:, 1+ i * pqNtor * MM : end])
        #display(r[:, 1: i * pqNtor * MM])
        
        #wot the actual fuck is going on here.
        r2D_alpha[1+idx * fM : (idx+1)*fM, 1 :(q-i) * N * MM] = r[:, 1+ i * N * MM : end]

        r2D_alpha[1+idx * fM : (idx+1)*fM, 1 + (q-i) * N * MM : end] = r[:, 1: i * N * MM]

        θ2D_alpha[1+idx * fM : (idx+1)*fM, 1 :(q-i) * N * MM] = θ[:, 1+i * N * MM : end]

        θ2D_alpha[1+idx * fM : (idx+1)*fM, 1 + (q-i) * N * MM : end] = θ[:, 1: i * N * MM]
    
    end

    #display(r2D_alpha) 
    #display(t2D_alpha)

    r2D_vartheta = zeros((qfM, MM * N))
    θ2D_vartheta = zeros((qfM, MM * N))

    for i in 0:MM * N-1
        #v odd that this is required. but otherwise the mod function removes some of the input???
        arr = 0:qfM
        idx = @. mod(arr - i * p, qfM) .+ 1

        


        idx = idx[1:end-1]
        #println(idx)

        #this will not work.
        r2D_vartheta[:, i+1] = r2D_alpha[idx, i+1]
        θ2D_vartheta[:, i+1] = θ2D_alpha[idx, i+1]
    end

    #display(t2D_vartheta)
    #odd to change r to s here but leave t as theta
    #this should be vartheta I think,
    #unclear what is going on though.
    scos_surf, ssin_surf = rfft2D(r2D_vartheta, M, N)
    tcos_surf, tsin_surf = rfft2D(θ2D_vartheta, M, N)

    return scos_surf, tsin_surf, ssin_surf, tcos_surf
    #ideally we will swap this order
    #it is pretty fkn stupid.
    #return scos_surf, θsin_surf, ssin_surf, θcos_surf

    #we can work on the other shite later.
end


    
#this is the equation of δS. which is a function of the fourier coeffiecients for r and θ. However, these are combined into the x variable for root finding.
#if we can fix the fft, this is probably as optimal as it will get tbh.
#the fft is a disaster, as we do it a lot and each time it is allocating a fair big of memoyr.
#function action_grad!(δS::Array{Float64}, x::Array{Float64}, p::It64, q::Int64, a::Float64, qN::Int64, ζ::Float64, nlist::Array{Int64}, prob::ProblemT, met::MetT, B::BFieldT, Br::Array{Float64}, Bθ::Array{Float64}, Bζ::Array{Float64})
#need some structs, this is a disaster.
#something about this is completly cooked
#cannot work out what.
#Beyond confused by this function, times to abandon.
function action_grad!(δS::Array{Float64}, x::Array{Float64}, CT::CoefficientsT, p::Int64, q::Int64, a::Float64, Ntor::Int64, r::Array{Float64}, θ::Array{Float64}, ζ::StepRangeLen{Float64}, ift_1D_p::AbstractFFTs.ScaledPlan, prob::ProblemT, met::MetT, B::BFieldT, Br::Array{Float64}, Bθ::Array{Float64}, Bζ::Array{Float64}, ift_r1D::Array{ComplexF64}, ift_θ1D::Array{ComplexF64}, nlist,rdot_fft_cos::Array{Float64}, rdot_fft_sin::Array{Float64}, θdot_fft_cos::Array{Float64}, θdot_fft_sin::Array{Float64}, ft_r1D::Array{ComplexF64}, ft_θ1D::Array{ComplexF64}, rft_1D_p::FFTW.rFFTWPlan)
#function action_grad!(δS::Array{Float64}, x::Array{Float64}, p::Int64, q::Int64, a::Float64, qN::Int64, ζ, nlist, prob::ProblemT)

    iota = p/q


    unpack_coeffs!(x, CT, Ntor)

    area = CT.θcos[1]
    #area = tcos[1]

    #r = irfft1D(rcos, rsin)
    

    #θ = irfft1D(tcos, tsin)

    get_r_t!(r, θ, CT, ift_1D_p, ift_r1D, ift_θ1D, Ntor)
    #get_r_t!(r, θ, rcos, tsin, rsin, tcos, ift_1D_p, ft_r1D, ft_θ1D, Ntor)
    

    #do not have this here.
    θ .+= iota .* ζ

   
    #this is unfortunatly inevitable.
    #think we will need to recompute the B field for this and for the gradient. Don't think it can be shared unfort.
    for i in 1:1:length(r)
        #perhaps r, t, z are not what we think they are? 
        #pretty hekin likely I think
        #nah cause the test compute B uses them as is?
        #these two functions are the main problemo now.
        #perhaps it is worth implementing a compute covarient B or something?
        #in this case we don't need all the extra shibang.
        #probably is worth it.
        #also, we will just need the Jacobian
        #this will also be a useful tool if we want to add poincare plots to our module.
        #we can probably just have a compute B where instead of passing in a BFieldT and a MetT we just pass in a 3 array and a single jacobian value.
        #could also be a problem that our structs are mutable
        #maybe worth changing met.J to met.J[1] etc
        prob.compute_met(met, r[i], θ[i], ζ[i], prob.geo.R0)

        #unsure why it is still unhappy.
        #maybe this is just how the q-profile works
        compute_B!(B.B, met.J, prob.q, prob.isl, prob.isl2, r[i], θ[i], ζ[i])

        #unsure how jacobian sits with all of this.
        Br[i] = B.B[1]
        Bθ[i] = B.B[2]
        Bζ[i] = B.B[3]
    end

    #to be fourier transformed.
    #these are still be allocated here!
    θdot = @. Bθ / Bζ

    
    #nv may be the biggest mystery atm.
    #it may be the lagrange multiplier that is mentioned?
    rdot = @. Br / Bζ - CT.nv[1] / Bζ

    
    #display(rhs_rdot)
    #println(rhs_rdot_fft_cos)


    rfft1D!(rdot_fft_cos, rdot_fft_sin, rdot, ft_r1D, rft_1D_p)
    rfft1D!(θdot_fft_cos, θdot_fft_sin, θdot, ft_θ1D, rft_1D_p)

    #rdot_fft_cos, rdot_fft_sin = rfft1D(rdot)

    #θdot_fft_cos, θdot_fft_sin = rfft1D(θdot)
    
    #println(rhs_tdot_fft_cos)
    
    #unfort that we are doing this.
    #ideally it would not be done.
    #at least not twice... (inside of get_r_θ)
    #rcos = view(x, 1:Ntor)
    #rsin = view(x, Ntor+1:2*Ntor)
    #tcos = view(x, 2*Ntor+1:3*Ntor)
    #tsin = view(x, 3*Ntor+1:4*Ntor)

    #display(tcos)

    #δS = zeros(4*Ntor + 1)
    
    #δS = [rcos ; rsin ; tcos ; rsin]

    #so we are dropping one from θcos and from rsin? seems odd.

    #ok lets revert the shape back

    

    #need to actually understand these equations.
    #δS[1:Ntor] = @. CT.rsin * nlist / q - rhs_rdot_fft_cos[1:Ntor]
    #δS[Ntor+1:2*Ntor] = @. (CT.θsin * nlist / q - rhs_tdot_fft_cos[1:Ntor])
    #δS[2*Ntor+1:3*Ntor-1] = @. (-CT.θcos * nlist / q - rhs_tdot_fft_sin[1:Ntor])[2:end]
    #δS[3*Ntor+1-1:4*Ntor-2] = @. (-CT.rcos * nlist / q - rhs_rdot_fft_sin[1:Ntor])[2:end]
    #δS[3*Ntor+1-1] += iota #still got no idea.
    #δS[end] = area - a #this is the lagrange multiplier constraint.
    #display(length(δS))
    #display(length(tcos))
    #qN = Ntor - 1
    #δS[1] = area - a

    
    #new packing
    #[rcos ; rsin[2:end] ; tcos ; tsin[2:end] ; [nv]]

    #wow we finally got this working 
    #holy moly.
    δS[1:Ntor] = @. CT.rsin * nlist / q - rdot_fft_cos[1: Ntor]


    δS[Ntor+1:2*Ntor] = @. (CT.θsin * nlist / q - θdot_fft_cos[1 : Ntor])

    δS[Ntor+1] += iota

    δS[2*Ntor+1:3*Ntor-1] = @. (-CT.θcos * nlist / q - θdot_fft_sin[1:Ntor])[2:end]

    δS[3*Ntor+1-1:4*Ntor-2] = @. (-CT.rcos * nlist / q - rdot_fft_sin[1:Ntor])[2:end]

    δS[end] = area - a

    #δS[Ntor+1:2*Ntor] = @. (CT.θsin * nlist / q - rhs_tdot_fft_cos[1 : qN + 1])

    #δS[Ntor+1] += iota

    #δS[2*Ntor+1:3*Ntor-1] = @. (-CT.θcos * nlist / q - rhs_tdot_fft_sin[1:qN+1])[2:end]

    #δS[3*Ntor+1-1:4*Ntor-1-1] = @. (-CT.rcos * nlist / q - rhs_rdot_fft_sin[1:qN+1])[2:end]

    #δS[end] = area - a
    
    #og packing
    #[[nv]; rcos ; tsin[2:end] ; rsin[2:end] ; tcos]
    #qN = Ntor - 1
    #δS[1] = area - a

    #δS[2 : qN + 2] = @. CT.rsin * nlist / q - rhs_rdot_fft_cos[1: qN + 1]

    #δS[qN+3: 2*qN + 2] = @. (-CT.rcos * nlist / q - rhs_rdot_fft_sin[1:qN+1])[2:end]

    #δS[2*qN+3 : 3*qN + 3] = @. (CT.θsin * nlist / q - rhs_tdot_fft_cos[1 : qN + 1])

    #δS[2 * qN + 3] += iota #wot.

    #δS[3 * qN + 4 : end] = @. (-CT.θcos * nlist / q - rhs_tdot_fft_sin[1:qN+1])[2:end]
    

end


function rfft1D!(cos_coef::Array{Float64}, sin_coef::Array{Float64}, func::Array{Float64}, fft_res::Array{ComplexF64}, plan::FFTW.rFFTWPlan)

    #shouldn't be doing this.
    Nfft = length(func)

    #test = rfft(func)
    #display(length(test))
    #display(length(func))
    #display(length(fft_res))

    #display(length(plan * func))

    #fft_res .= rfft(func)

    fft_res .= plan * func

    cos_coef .= real.(fft_res) ./ Nfft .* 2

    cos_coef[1] /= 2

    sin_coef .= @. -imag(fft_res) / Nfft * 2

    sin_coef[1] = 0.0


end

#reconstructs r and t from the combined fourier coefficients.
#ideally this will have a a plan in here as well.
#takes in inverse real fourier transform to get r, θ from the coefficients.
function get_r_t!(r::Array{Float64}, θ::Array{Float64}, CT::CoefficientsT, ft_plan::AbstractFFTs.ScaledPlan, ft_r1D::Array{ComplexF64}, ft_θ1D::Array{ComplexF64}, Ntor::Int64)

    #probably need to know the size of this already!
    #nfft_multiplier = 2
    #Nfft = nfft_multiplier * (size(cos_in)[end] - 1)

    #2 is an optional arg.
    Nfft = 2 * (Ntor - 1)

    #lets change the structure of x.
    #each has dimension qN + 1
    #no removing the fkn zeros or anything stupid anymore.
    #think we want to work with Npol and Ntor as our 2 actual dimensions
    #Ntor = qN + 1
    #rcos = view(x, 1:Ntor)
    #rsin = view(x, Ntor+1:2*Ntor)
    #tcos = view(x, 2*Ntor+1:3*Ntor)
    #tsin = view(x, 3*Ntor+1:4*Ntor)

    #display(rcos)

    #ensures how fft array is of a minimum length, ideally this would be done.
    #I think ft_pad could be worked out in advance, then 
    #we just have the extra values always set to zero?
    #should always be the case.
    #Npadr = floor(Int64,  Nfft / 2 - length(rcos)/2) + 1

    #Npadt = floor(Int64,  Nfft / 2 - length(tcos)/2) + 1
    #mayhap we don't actually need ft_padr, as if we do it in place
    #we can just use r?
    ft_r1D[1:Ntor] .= @. (CT.rcos - 1im * CT.rsin) * Nfft
    ft_θ1D[1:Ntor] .= @. (CT.θcos - 1im * CT.θsin) * Nfft

    #display(rcos)
    #display(Nfft)
    #display(ft_r1D)

    ft_r1D[1] *= 2
    ft_θ1D[1] *= 2
    #unsure exactly what is going on here tbh.
    r .= ft_plan * ft_r1D
    θ .= ft_plan * ft_θ1D


    
end


function get_r_t!(r::Array{Float64}, θ::Array{Float64}, rcos, θsin, rsin, θcos, ft_plan::AbstractFFTs.ScaledPlan, ft_r1D::Array{ComplexF64}, ft_θ1D::Array{ComplexF64}, Ntor::Int64)

    #probably need to know the size of this already!
    #nfft_multiplier = 2
    #Nfft = nfft_multiplier * (size(cos_in)[end] - 1)

    #2 is an optional arg.
    Nfft = 2 * (Ntor - 1)

    #lets change the structure of x.
    #each has dimension qN + 1
    #no removing the fkn zeros or anything stupid anymore.
    #think we want to work with Npol and Ntor as our 2 actual dimensions
    #Ntor = qN + 1
    #rcos = view(x, 1:Ntor)
    #rsin = view(x, Ntor+1:2*Ntor)
    #tcos = view(x, 2*Ntor+1:3*Ntor)
    #tsin = view(x, 3*Ntor+1:4*Ntor)

    #display(rcos)

    #ensures how fft array is of a minimum length, ideally this would be done.
    #I think ft_pad could be worked out in advance, then 
    #we just have the extra values always set to zero?
    #should always be the case.
    #Npadr = floor(Int64,  Nfft / 2 - length(rcos)/2) + 1

    #Npadt = floor(Int64,  Nfft / 2 - length(tcos)/2) + 1
    #mayhap we don't actually need ft_padr, as if we do it in place
    #we can just use r?
    ft_r1D[1:Ntor] .= @. (rcos - 1im * rsin) * Nfft
    ft_θ1D[1:Ntor] .= @. (θcos - 1im * θsin) * Nfft

    #display(rcos)
    #display(Nfft)
    #display(ft_r1D)

    ft_r1D[1] *= 2
    ft_θ1D[1] *= 2
    #unsure exactly what is going on here tbh.
    r .= ft_plan * ft_r1D
    θ .= ft_plan * ft_θ1D


    
end




function pack_coeffs!(x::Array{Float64}, CT::CoefficientsT, N)

    #I think we will still want to change the order of this.
    #just need to be a wee bit careful.
    #should be easier with this function and the struct.
    #x[1] = CT.nv[1]
    #changed the order hopefully not a mistake!
    #x[1:N] = CT.rcos[:]
    #x[N+1:2*N-1] = CT.rsin[2:end]
    #x[2*N-1+1:3*N-1] = CT.θcos[:]
    #x[3*N-1+1:4*N-1-1] = CT.θsin[2:end]
    #x[end] = CT.nv[1]
    #so the plus 1 is extremely important for efficient solving!
    #why??>>?>
    #og packing
    #x .= [CT.nv ; CT.rcos ; CT.θsin[2:end] ; CT.rsin[2:end] ; CT.θcos] .+ 1
    #new packing
    x .= [CT.rcos ; CT.rsin[2:end] ; CT.θcos ; CT.θsin[2:end] ; CT.nv] .+ 1

end

function unpack_coeffs!(x::Array{Float64}, CT::CoefficientsT, N)

    #I think we will still want to change the order of this.
    #just need to be a wee bit careful.
    #should be easier with this function and the struct.
    #x[1] = CT.nv[1]
    
    #changed the order hopefully not a mistake!
    CT.rcos[:] = x[1:N] .- 1
    CT.rsin[2:end] = x[N+1:2*N-1] .- 1
    CT.θcos[:] = x[2*N-1+1:3*N-1] .-1
    CT.θsin[2:end] = x[3*N-1+1:4*N-1-1] .-1
    CT.nv[1] = x[end] - 1
    CT.rsin[1] = 0.0
    CT.θsin[1] = 0.0

    #og packing
    #=
    qN = N - 1
    CT.nv[1] = x[1] - 1
    CT.rcos .= x[2:qN + 2]  .- 1

    CT.θsin .= [[0] ; x[qN + 3:2 * qN + 2] .- 1]
    CT.rsin .= [[0] ; x[2*qN + 3 : 3 * qN + 2] .- 1]
    CT.θcos .= x[3 * qN + 3 : end] .- 1
    =#

    

end


#computes the jacobian matrix for the action grad for more efficient root finding!
function action_grad_jm()

end