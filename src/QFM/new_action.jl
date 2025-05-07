#combining the two is almost certainly not worth the hassle.
#this may not be the best use of time.
struct NewCoefficientsT
    nv :: Array{Float64, 1} #this is actually just a float. but we want to mutate it.
    rcos :: Array{Float64, 1}
    θcos :: Array{Float64, 1}
    rsin :: Array{Float64, 1}
    θsin :: Array{Float64, 1}
    function CoefficientsT(N::Int64)
        new(zeros(1), zeros(N), zeros(N), zeros(N), zeros(N))
    end
end

#this file is probably about as optimized as we can hope for,
#we could perhaps add a fj! function to combine, 
#but this file needs a serios clearn it is completly fked.

#ok, to summarise the fourier transform sizes
#our actual data, the coeffs are length qN  
#we expand this with a fft_alias, so our arrays are
#nfft * qN + 1  (real's)(+1 is the n=0 term, only needed for cos terms)
#the coeff arrays are irfft'ed
#irfft arrays have size 2*(N-1), where N is the size of the rfft.
#this gives r, θ
#which must be (nfft*qN) * 2 in length (+1 and -1 cancel out)
#this is the size of most of our working data, eg B etc.
#therefore, most of our data will either be 
#nfft*qN + 1
#or
#2*nfft*qN
#despite the actual coeff arrays being qN + 1 in size.


#struct to condense the inputs for action grad. allowing preallocation of variables.
#perhaps this should be split up into multiple structs. this is wild.
#the split between this struct and the next is much more arbitrary than we would like.
struct NewActionGradInputsT
    a :: Int64
    b :: Int64
    qN :: Int64
    nfft :: Int64
    r :: Array{Float64}
    θ :: Array{Float64}
    ζ :: StepRangeLen{Float64}
    prob :: ProblemT
    met :: MetT
    B :: BFieldT
    nlist :: Array{Int64}
    rdot_rhs :: Array{Float64}
    θdot_rhs :: Array{Float64}
    rdot_rhs_cos :: Array{Float64} #shouldn't be here
    rdot_rhs_sin :: Array{Float64}
    θdot_rhs_cos :: Array{Float64}
    θdot_rhs_sin :: Array{Float64}
    function NewActionGradInputsT(a::Int64, b::Int64, qN::Int64, nfft::Int64, ζ::StepRangeLen{Float64}, prob::ProblemT, met::MetT, B::BFieldT, nlist::UnitRange{Int64})
        #size of the ifft'd arrays.
        arr_size = 2 * qN * nfft
        irfft_size = qN*nfft+1
        new(a, b, qN, nfft, zeros(arr_size), zeros(arr_size), ζ, prob, met, B, nlist, zeros(arr_size), zeros(arr_size), zeros(irfft_size), zeros(irfft_size), zeros(irfft_size), zeros(irfft_size))
    end
end

#ideally, this will mainly just be plans, perhaps a temp array or 2.
struct NewFTDataT
    irfft_1D_p :: AbstractFFTs.ScaledPlan
    rfft_1D_p :: FFTW.rFFTWPlan
    temp :: Array{ComplexF64}
    jm_temp :: Array{ComplexF64, 2} #probably 2 dims.
    jmrfft_1D_p :: FFTW.rFFTWPlan #will be tricky af.
    function NewFTDataT(qN::Int64, nfft::Int64)
        #size of arrays in fourier space
        fft_size = 2 * qN * nfft
        temp = zeros(fft_size)
        rfft_1D_p = plan_rfft(temp)
        #size of arrays after they have been irfft'd.
        ifft_size = qN * nfft + 1
        temp = zeros(ComplexF64, ifft_size)
        irfft_1D_p = plan_irfft(temp, fft_size)

        comb_size = 4 * qN + 2
        temp = zeros(fft_size, comb_size)
        jmrfft_1D_p = plan_rfft(temp, [1]) #this may need to change!
        new(irfft_1D_p, rfft_1D_p, zeros(ifft_size), zeros(ifft_size, comb_size), jmrfft_1D_p)
        #new(irft_1D_p, rft_1D_p, zeros(ComplexF64, fft_size+1), zeros(ComplexF64, fft_size+1), zeros(fft_size+1), zeros(fft_size+1), zeros(fft_size+1), zeros(fft_size+1), zeros(ComplexF64, fft_size+1), zeros(ComplexF64, fft_size+1))
    end
end


#struct to condense the inputs for action grad jm. allowing preallocation of variables.
#this takes exactly the same args. Don't think it should!
#still needs work, we will need to fix action_grad_jm first tbh.
#this is even more of a disaster. Probably want the normal one and an extra one, just so this isn't quite so fkn cooked.
struct NewActionGradInputsJMT
    ν_term :: Array{Float64}
    ν_cos :: Array{Float64}
    ν_sin :: Array{Float64}
    cosnzq :: Array{Float64, 2}
    sinnzq :: Array{Float64, 2}
    rdot_rhs_drc :: Array{Float64, 2}
    rdot_rhs_drs :: Array{Float64, 2}
    rdot_rhs_dθc :: Array{Float64, 2}
    rdot_rhs_dθs :: Array{Float64, 2}
    θdot_rhs_drc :: Array{Float64, 2}
    θdot_rhs_drs :: Array{Float64, 2}
    θdot_rhs_dθc :: Array{Float64, 2}
    θdot_rhs_dθs :: Array{Float64, 2}
    drdot :: Array{Float64, 2}
    dθdot :: Array{Float64, 2}
    drdot_cos ::Array{Float64, 2}
    drdot_sin ::Array{Float64, 2}
    dθdot_cos ::Array{Float64, 2}
    dθdot_sin ::Array{Float64, 2}
    function NewActionGradInputsJMT(ip::NewActionGradInputsT)
        arr_size = 2 * ip.qN * ip.nfft
        comb_size = 4*ip.qN+2
        irfft_size = ip.qN*ip.nfft+1
        nzq = ip.ζ * ip.nlist' / ip.a
        #this has got to be wrong lol.
        new(zeros(arr_size), zeros(irfft_size), zeros(irfft_size), cos.(nzq), sin.(nzq), zeros(arr_size, ip.qN+1), zeros(arr_size, ip.qN+1), zeros(arr_size, ip.qN+1), zeros(arr_size, ip.qN+1), zeros(arr_size, ip.qN+1), zeros(arr_size, ip.qN+1), zeros(arr_size, ip.qN+1), zeros(arr_size, ip.qN+1), zeros(arr_size, comb_size), zeros(arr_size, comb_size), zeros(irfft_size, comb_size), zeros(irfft_size, comb_size), zeros(irfft_size, comb_size), zeros(irfft_size, comb_size))
    end
end



function new_action(a::Int64, b::Int64, prob::ProblemT, met::MetT, B::BFieldT, M::Int64, N::Int64, sguess=0.5::Float64, nfft=2::Int64)
    #a, b is the toroidal, poloidal mode number that defines the rational surfaces a/b
    #surfies is defined at q=a/b
    #not that b≡p and a≡q in original case, as they use iota.

    #N is the number of fourier harmonics included in the trial function.
    #this gets expanded to q*N, as the new field line will have 2πq periodicity
    #so the domain is expanded, requiring more fourier harmonics.
    #placeholders until we swap fully.
    q = a
    p = b
    qN = q * N

    MM = 2 * nfft #still very unsure what this is.

    #number of pseudo field lines to be found
    fl = MM * N

    qfM = fl * q
    Nfft = qfM
    #r(ζ) = r_0^c + ∑_{n=1}^qN (r_n^c cos(nζ/q) + r_n^s sin(nζ/q))
    #θ(ζ) = θ_0^c + ι*ζ + ∑_{n=1}^qN (θ_n^c cos(nζ/q) + θ_n^s sin(nζ/q))


    dζ = 2π / fl
    dθ = 2π / qfM

    nlist = range(0, qN)
    #don't think zeta is the toroidal component.
    ζ = range(0, Nfft-1) .* dζ

    #array storing the unknowns
    #[[ν], [rcos], [θsin[2:end]], [rsin[2:end]], [θcos]]
    #sin[1] components are ignored as sin(0) = 0.
    #qN counts from n=1 to n=qN
    #4 * qN is the 4 set of coeffs for n=1 to n=qN
    #+1 for ν
    #+2 for the two cos terms which are n=0

    x0 = zeros(4*qN + 1 + 2)

    #these store the minimised values for the coefficeints at each iteration
    rcosarr = zeros(fl , qN+1)
    rsinarr = zeros(fl , qN+1) #does not actually need th +1 (n=0) but kept for symmetry.
    θcosarr = zeros(fl , qN+1)
    θsinarr = zeros(fl , qN+1)
    nvarr = zeros(fl)

    #struct storing the inputs needed for action grad
    ip = NewActionGradInputsT(a, b, qN, nfft, ζ, prob, met, B, nlist)

    #creates struct for storing fft plans and temp arrays
    #for efficient fourier transforming
    ftd = NewFTDataT(qN, nfft)

    #this will be so wrong lol
    ipjm = NewActionGradInputsJMT(ip)


    coefs = CoefficientsT(qN+1)

    #iterates over different poloidal areas, this creates multiple psuedo fieldlines that minamise the action.
    #these are then combined into a surface below.
    for jpq in 0:fl-1
        #maybe this is supposed to be alpha?
        #or area? fk knows.
        a = jpq * dθ
        if jpq == 0
            #initial guesses for the alg.
            coefs.rcos[1] = sguess
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
            coefs.rcos .= rcosarr[jpq, :]
            coefs.rsin .= rsinarr[jpq, :]
            coefs.θcos .= θcosarr[jpq, :]
            coefs.θsin .= θsinarr[jpq, :]
            coefs.nv[1] = nvarr[jpq]
            rcos0 = rcosarr[jpq, :]
            rsin0 = rsinarr[jpq, :]
            θcos0 = θcosarr[jpq, :]
            θsin0 = θsinarr[jpq, :]
            nv0 = nvarr[jpq]

            θcos0[1]  += dθ 
            coefs.θcos[1] += dθ
        end
        #combine the values into a single array for the root findning alg's argument
        new_pack_coeffs!(x0, coefs)

        #creates a standin function that fits the form needed for NLSolve.
        ag!(δS, x) = new_action_grad!(δS, x, a, coefs, ip, ftd)

        #creates a standin function that fits the form needed for NLSolve.
        #this function defines the gradient, improving efficiency.
        ag_JM!(JM, x) = new_action_grad_jm!(JM, x, a, coefs, ip, ftd, ipjm)

        #solve for the psuedo field line.
        sol = nlsolve(ag!, ag_JM!, x0)
        #sol = nlsolve(ag!, x0)

        #unpacks the single array solution into CoefficientsT struct.
        new_unpack_coeffs!(sol.zero, coefs, qN+1)

        #stores this field line.
        rcosarr[jpq+1, :] = coefs.rcos
        rsinarr[jpq+1, :] = coefs.rsin
        θcosarr[jpq+1, :] = coefs.θcos
        θsinarr[jpq+1, :] = coefs.θsin
        nvarr[jpq+1] = coefs.nv[1]
        
    end

    return wrap_field_lines(rcosarr, rsinarr, θcosarr, θsinarr, MM, M, N, p, q, qfM, Nfft)
end

#this one computes the action gradient and the second derivative.
#this will also hopefully be a bit clearer.
function new_action_grad!(δS::Array{Float64}, x::Array{Float64}, a::Float64, coefs::CoefficientsT, ip::NewActionGradInputsT, ftd::NewFTDataT)
   
    #change to be just in terms of qN
    new_unpack_coeffs!(x, coefs, ip.qN+1)

    #should be changed to b / a to match other work.
    iota = ip.b / ip.a

    #note this is the area of the trial function
    #not sure it is actually clearer to have this as its own variable.
    area = coefs.θcos[1]


    #inverse fourier transform, done for r and θ at the same time.
    #problemo is here!
    irfft1D!(ip.r, coefs.rcos, coefs.rsin, ftd.irfft_1D_p, ftd.temp, ip.qN, ip.nfft)
    irfft1D!(ip.θ, coefs.θcos, coefs.θsin, ftd.irfft_1D_p, ftd.temp, ip.qN, ip.nfft)
    #get_r_t!(ip.r, ip.θ, coefs, ftd.ift_1D_p, ftd.ift_r1D, ftd.ift_θ1D, ip.Ntor)

    #ip.r .= irfft1D(coefs.rcos, coefs.rsin)
    #ip.θ .= irfft1D(coefs.θcos, coefs.θsin)

    #reflects that there is an additional term in the trial function for θ
    ip.θ .+= iota .* ip.ζ


    for i in 1:2*ip.qN*ip.nfft
        ip.prob.met(ip.met, ip.r[i], ip.θ[i], ip.ζ[i], ip.prob.geo.R0)

        compute_B!(ip.B, ip.met, ip.prob.q, ip.prob.isls, ip.r[i], ip.θ[i], ip.ζ[i])

        #these are awful names,
        #θdot is what θdot is equal to, when the gradient is zero.
        #this is why Zhisong used rhs terminology.
        #these are essentially the parts of the two equations that need to be fourier transformed.
        ip.θdot_rhs[i] = ip.B.B[2] / ip.B.B[3]

        #this is fine
        #probably should have a /2π q in ν term for full consistency.
        #to be certain on this, we will need to understand the scaling for fft,
        #i.e. does that match the analytical case of integration.
        #if the action starts with -ν this should be -ν here.
        ip.rdot_rhs[i] = ip.B.B[1] / ip.B.B[3] - coefs.nv[1] / (2*ip.a*π*ip.B.B[3] * ip.met.J[1])
        #ip.rdot_rhs[i] = ip.B.B[1] / ip.B.B[3] - coefs.nv[1] / (ip.B.B[3] * ip.met.J[1])
    end

    #this is fine
    rfft1D!(ip.rdot_rhs_cos, ip.rdot_rhs_sin, ip.rdot_rhs, ftd.temp, ftd.rfft_1D_p, 2*ip.nfft*ip.qN)
    rfft1D!(ip.θdot_rhs_cos, ip.θdot_rhs_sin, ip.θdot_rhs, ftd.temp, ftd.rfft_1D_p, 2*ip.nfft*ip.qN)




    #equations are structured so that first set comes from equation with rdot (i.e. δS/δθ) multiplied by cos
    #because cos(n=0) ≠ 0, we have aN + 1 equations.
    #second set comes from same equation multiplied by sin, giving aN equations
    #third set come from θdot eqaution (i.e. δS/δr) multiplied by cos
    #again, this has aN + 1 equations
    #fourth set is same equations with sin
    #giving aN equations
    #final set is a single equation for ν.
    #[rcos ; rsin[2:end] ; tcos ; tsin[2:end] ; [nv]]
    #δS has length 4*qN + 3.
    #[qN+1 ; qN ; qN+1 ; qN ; 1]
    #rcos coefs come from δS/δθ * cos
    δS[1:ip.qN+1] = @. coefs.rsin * ip.nlist / ip.a - ip.rdot_rhs_cos[1:ip.qN+1]

    #[2:end] as we dont need the rsin[1] coef, as it is zero
    #rsin coefs come from δS/δθ * sin
    #only qN coeffs for r^s
    δS[ip.qN+2:2*ip.qN+1] = @. (-coefs.rcos * ip.nlist / ip.a - ip.rdot_rhs_sin[1:ip.qN+1])[2:end]

    #θcos
    δS[2*ip.qN+2:3*ip.qN+2] = @. coefs.θsin * ip.nlist / ip.a - ip.θdot_rhs_cos[1:ip.qN+1]

    #(n=0) term for cos has additional iota term
    δS[2*ip.qN+2] += iota

    #θsin
    δS[3*ip.qN+3:4*ip.qN+2] = @. (-coefs.θcos * ip.nlist / ip.a - ip.θdot_rhs_sin[1:ip.qN+1])[2:end]

    #lagrange multiplier term
    δS[end] = area - a



    #note that we are essentially following Stuarts code for this, Zhisong's could be 
    #a bit wrong
    #δS[1] = area - a

    #equation defining the r^c coeffs
    #comes from δS/δθ * cos(nζ/q) equation
    #rsin comes from ṙ.
    #δS[2:qN+1] = @. coeffs.rsin * ip.nlist / ip.a - ftd.rdot_fft_cos

    #equation defining the θ^s coefs
    #comes from δS/δr * sin equation
    #r
    #δS[qN+3: 2*qN + 2] = @. (-coefs.rcos * ip.nlist / ip.q - ftd.rdot_fft_sin[1:qN+1])[2:end]

    #δS[2*qN+3 : 3*qN + 3] = @. (coefs.θsin * ip.nlist / ip.q - ftd.θdot_fft_cos[1 : qN + 1])

    #δS[2 * qN + 3] += iota #wot.

    #δS[3 * qN + 4 : end] = @. (-coefs.θcos * ip.nlist / ip.q - ftd.θdot_fft_sin[1:qN+1])[2:end]
end

function new_action_grad_jm!(JM::Array{Float64, 2}, x::Array{Float64}, a, coefs::CoefficientsT, ip::NewActionGradInputsT, ftd::NewFTDataT, ipjm::NewActionGradInputsJMT)

    new_unpack_coeffs!(x, coefs, ip.qN+1)

    #should be changed to b / a to match other work.
    iota = ip.b / ip.a

    #note this is the area of the trial function
    area = coefs.θcos[1]


    #inverse fourier transform, done for r and θ at the same time.
    irfft1D!(ip.r, coefs.rcos, coefs.rsin, ftd.irfft_1D_p, ftd.temp, ip.qN, ip.nfft)
    irfft1D!(ip.θ, coefs.θcos, coefs.θsin, ftd.irfft_1D_p, ftd.temp, ip.qN, ip.nfft)
    #get_r_t!(ip.r, ip.θ, coefs, ftd.ift_1D_p, ftd.ift_r1D, ftd.ift_θ1D, ip.Ntor)

    #ip.r .= irfft1D(coefs.rcos, coefs.rsin)
    #ip.θ .= irfft1D(coefs.θcos, coefs.θsin)

    #reflects that there is an additional term in the trial function for θ
    ip.θ .+= iota .* ip.ζ

    for i in 1:2*ip.qN*ip.nfft
        ip.prob.met(ip.met, ip.r[i], ip.θ[i], ip.ζ[i], ip.prob.geo.R0)

        #need B.dB now, so we just compute the full thing!
        compute_B!(ip.B, ip.met, ip.prob.q, ip.prob.isls, ip.r[i], ip.θ[i], ip.ζ[i])

        rdot_rhs_dr = (ip.B.dB[1, 1]/ip.B.B[3] - ip.B.B[1] * ip.B.dB[3, 1] / ip.B.B[3]^2 
                       + coefs.nv[1] * (ip.B.dB[3, 1] / (2*ip.a*π*ip.B.B[3]^2*ip.met.J[1]) 
                                        + ip.met.dJ[1] / (2*ip.a*π*ip.B.B[3] * ip.met.J[1]^2)))
        rdot_rhs_dθ = (ip.B.dB[1, 2]/ip.B.B[3] - ip.B.B[1] * ip.B.dB[3, 2] / ip.B.B[3]^2 
                       + coefs.nv[1] * (ip.B.dB[3, 2] / (2*ip.a*π*ip.B.B[3]^2*ip.met.J[1]) 
                                        + ip.met.dJ[2] / (2*ip.a*π*ip.B.B[3] * ip.met.J[1]^2)))

        θdot_rhs_dr = ip.B.dB[2, 1]/ip.B.B[3] - ip.B.B[2] * ip.B.dB[3, 1] / ip.B.B[3]^2

        θdot_rhs_dθ = ip.B.dB[2, 2]/ip.B.B[3] - ip.B.B[2] * ip.B.dB[3, 2] / ip.B.B[3]^2

        #terrible name, this was the ooBζ
        #weird shape is so this fits in array below
        ipjm.ν_term[i] = 1 / (2*ip.a*π * ip.B.B[3] * ip.met.J[1])
        #ipjm.ν_term[i, 1] = 1 / ip.B.B[3]
        
        #rdot_rhs_dr = ip.B.dB[1, 1] / ip.B.B[3] - (ip.B.B[1] - coefs.nv[1]) / ip.B.B[3] * ip.B.dB[3, 1]
        #dBrdθ/Bζ - (Br - nv) / Bζ^2 * dBζdθ
        #rdot_rhs_dθ = ip.B.dB[1, 2] / ip.B.B[3] - (ip.B.B[1] - coefs.nv[1]) / ip.B.B[3] * ip.B.dB[3, 2]
        #dBθdr/Bζ - Bθ / Bζ^2 * dBζdr
        #θdot_rhs_dr = ip.B.dB[2, 1] / ip.B.B[3] - ip.B.B[2] / ip.B.B[3] * ip.B.dB[3, 1]
        #dBθdθ/Bζ - Bθ / Bζ^2 * dBζdθ
        #θdot_rhs_dθ = ip.B.dB[2, 2] / ip.B.B[3] - ip.B.B[2] / ip.B.B[3] * ip.B.dB[3, 2]
        
        #ipjm.ooBζ[i] = 1 / ip.B.B[3]


        #iterates over the n's
        for j in 1:ip.qN+1
            #this notation is a bit nicer, probbaly change everything from rcos to rc.
            #although having rs could be confusing with our new s variable.
            #but this doesn't come up heaps tbh.
            #derivative w.r.t r^c
            ipjm.rdot_rhs_drc[i, j] = rdot_rhs_dr * ipjm.cosnzq[i, j]
            #derivative w.r.t r^s
            ipjm.rdot_rhs_drs[i, j] = rdot_rhs_dr * ipjm.sinnzq[i, j]
            #derivative w.r.t θ^c
            ipjm.rdot_rhs_dθc[i, j] = rdot_rhs_dθ * ipjm.cosnzq[i, j]
            #derivative w.r.t θ^s
            ipjm.rdot_rhs_dθs[i, j] = rdot_rhs_dθ * ipjm.sinnzq[i, j]

            #derivative w.r.t r^c
            ipjm.θdot_rhs_drc[i, j] = θdot_rhs_dr * ipjm.cosnzq[i, j]
            #derivative w.r.t r^s
            ipjm.θdot_rhs_drs[i, j] = θdot_rhs_dr * ipjm.sinnzq[i, j]
            #derivative w.r.t θ^c
            ipjm.θdot_rhs_dθc[i, j] = θdot_rhs_dθ * ipjm.cosnzq[i, j]
            #derivative w.r.t θ^s
            ipjm.θdot_rhs_dθs[i, j] = θdot_rhs_dθ * ipjm.sinnzq[i, j]
        end
    end

    #here is where we would pack them all together like an insance person.
    #perhaps we will not do that?
    #so ideall, drdot, will just be in the the loop above, saving allocations.o
    #for now,
    #so this is deliberatly packed the same way, then we can efficeintly fourier trasnfrom it all
    #then it can be combined into the main array v easily.
    #drdot[end] is the dummy nv, we should be able to set that to zero right from the start
    #then nothing will ever change it.
    #ipjm.drdot .= [ipjm.rdot_rhs_drc  ipjm.rdot_rhs_drs[:, 2:end]  ipjm.rdot_rhs_dθc  ipjm.rdot_rhs_dθs[:, 2:end]  ipjm.ν_term]
    #ipjm.dθdot .= [ipjm.θdot_rhs_drc  ipjm.θdot_rhs_drs[:, 2:end]  ipjm.θdot_rhs_dθc  ipjm.θdot_rhs_dθs[:, 2:end]  ipjm.ν_term]
    #this fixed it, so our fourier transform of the extra bit isn't looking good, hard to know why though tbh.
    #could be because the 1st terms are being set to zero etc? unsure how to fix.
    #the fourier transform could just be wrong for the 1d array tbh.
    #note that this version does seem to be significantly quicker.
    #dum_ν = zeros(2*ip.qN*ip.nfft, 1)
    #ipjm.drdot .= [ipjm.rdot_rhs_drc  ipjm.rdot_rhs_drs[:, 2:end]  ipjm.rdot_rhs_dθc  ipjm.rdot_rhs_dθs[:, 2:end]  dum_ν]
    #ipjm.dθdot .= [ipjm.θdot_rhs_drc  ipjm.θdot_rhs_drs[:, 2:end]  ipjm.θdot_rhs_dθc  ipjm.θdot_rhs_dθs[:, 2:end]  dum_ν]
    ipjm.drdot .= [ipjm.rdot_rhs_drc  ipjm.rdot_rhs_drs[:, 2:end]  ipjm.rdot_rhs_dθc  ipjm.rdot_rhs_dθs[:, 2:end]]
    ipjm.dθdot .= [ipjm.θdot_rhs_drc  ipjm.θdot_rhs_drs[:, 2:end]  ipjm.θdot_rhs_dθc  ipjm.θdot_rhs_dθs[:, 2:end]]

    #take the fourier transform along the second dimension, to efficiently fourier trasnform everthing at once.
    #TODO
    rfft1D!(ipjm.drdot_cos, ipjm.drdot_sin, ipjm.drdot, ftd.jm_temp, ftd.jmrfft_1D_p, 2*ip.nfft*ip.qN)
    rfft1D!(ipjm.dθdot_cos, ipjm.dθdot_sin, ipjm.dθdot, ftd.jm_temp, ftd.jmrfft_1D_p, 2*ip.nfft*ip.qN)

    #dummy_ν = zeros(1, 2*ip.qN*ip.nfft)

    
    #drdot = [transpose(ipjm.rdot_rhs_drc) ;  transpose(ipjm.rdot_rhs_drs[:, 2:end]) ;  transpose(ipjm.rdot_rhs_dθc) ;  transpose(ipjm.rdot_rhs_dθs[:, 2:end]) ;  dummy_ν]
    #dθdot = [transpose(ipjm.θdot_rhs_drc) ; transpose(ipjm.θdot_rhs_drs[:, 2:end]) ; transpose(ipjm.θdot_rhs_dθc) ; transpose(ipjm.θdot_rhs_dθs[:, 2:end]) ; dummy_ν]

    #drdot_cos, drdot_sin = rfft1D_JM(drdot)
    #dθdot_cos, dθdot_sin = rfft1D_JM(dθdot)
    #fourier transform along the second dimension, gives us what we want, all at hte same time.

    #ipjm.drdot_cos .= drdot_cos
    #ipjm.drdot_sin .= drdot_sin
    #ipjm.dθdot_cos .= dθdot_cos
    #ipjm.dθdot_sin .= dθdot_sin

    #does the ν term derivatives.
    #ideally, this would all be put into the grand drdot, but that doesn't seem to work.
    rfft1D!(ipjm.ν_cos, ipjm.ν_sin, ipjm.ν_term, ftd.temp, ftd.rfft_1D_p, 2*ip.nfft*ip.qN)

    #ooBcos, ooBsin = rfft1D_simple(ipjm.ν_term[:, 1])

    #not sure if we do actually need to do this, most things should get set.
    #JM .= 0.0

    #drdot_cos, drdot_sin = rfft1D_JM(ipjm.drdot)

    #we have chosen to organise our variables such that
    #[rcos ; rsin[2:end] ; tcos ; tsin[2:end] ; [nv]]
    #so even though our equations can be arbitrarly organised, the derivatives cannot
    #i.e. JM[1, 1] must be the derivative of whatever the first equation was, with respect to our first variable, r_0^c.
    #etc
    #however, the choice of equation order doesn't actually matter.
    #as the equation don't correspond to any variable, they just appear to, becuase each equation only tneds to contain one.
    
    JM .= 0.0

    #[rcos ; rsin[2:end] ; tcos ; tsin[2:end] ; [nv]]

    #Derivtive of first set of equations, w.r.t every other coefficient.
    #end - 1 here as the new drdot doesn't include the ν terms.
    JM[1:ip.qN+1, 1:end-1] .= -ipjm.drdot_cos[1:ip.qN+1, 1:end]
    #JM[1:ip.qN+1, 1:end] .= -drdot_cos[1:ip.qN+1, 1:end]

    #display(length(ip.nlist))
    #display(ip.qN)
    #display(length(CartesianIndex.(1:ip.qN+1, ip.qN+2:2*ip.qN+1)))
    #we are excluding the n=0 contribution
    JM[CartesianIndex.(2:ip.qN+1, ip.qN+2:2*ip.qN+1)] += ip.nlist[2:end] / ip.a

    #same with r^s terms, note that we still don't include the sin(0) term.
    JM[ip.qN+2:2*ip.qN+1, 1:end-1] .= -ipjm.drdot_sin[2:ip.qN+1, 1:end]

    JM[CartesianIndex.(ip.qN+2:2*ip.qN+1, 2:ip.qN+1)] += -ip.nlist[2:end]/ip.a

    #θ^c terms
    JM[2*ip.qN+2:3*ip.qN+2, 1:end-1] .= -ipjm.dθdot_cos[1:ip.qN+1, 1:end]

    JM[CartesianIndex.(2*ip.qN+3:3*ip.qN+2, 3*ip.qN+3:4*ip.qN+2)] += ip.nlist[2:end]/ip.a

    #θ^s terms
    JM[3*ip.qN+3:4*ip.qN+2, 1:end-1] .= -ipjm.dθdot_sin[2:ip.qN+1, 1:end]

    #fourth set with respect to θ^c (no θ^s derivs)
    JM[CartesianIndex.(3*ip.qN+3:4*ip.qN+2, 2*ip.qN+3:3*ip.qN+2)] += -ip.nlist[2:end]/ip.a

    #extra term from ν equation
    #this is d/dθ_0^c.
    JM[end, 2*ip.qN+2] += 1.0

    #in theory, these are not needed anymore, as there are in the big fourier transform thingo.
    JM[1:ip.qN+1, end] += ipjm.ν_cos[1:ip.qN+1]
    JM[ip.qN+2:2*ip.qN+1, end] += ipjm.ν_sin[2:ip.qN+1]


    #first set of eqns
    #[1:qN+1]
    #second set (only qN)
    #[qN+2:2*qN+1]
    #third set
    #[2*qN+2:3*qN+2]
    #fourth set (only qN)
    #[3*qN+3:4*qN+2]
    #fifth set (single equation for ν)
    #[4*qN+3] or just [end]

    #now we add the additional terms that come from the derivatives 
    #occuring in the expansion for r and θ
    #these derivatives pick out a specific trig term in the expansion
    #which integrate to 1, however, we have the n/q term from the derivative w.r.t ζ
    #i.e. the rdot and θdot terms
    #note that these indicies are pairwise, as the derivative are with respect to a specific n
    #so only a single equation is included

    #TODO Need to fix the signs here and probably everywhere tbh.
    #these indexes also won't work, as we need to avoid the sin_0 terms
    #first set with respect to r^s (no r^c derivs)
    #JM[CartesianIndex.(1:ip.qN+1, ip.qN+2:2*ip.qN+1)] += ip.nlist/ip.a
    #second set with respect to r^c (no r^s derivs)
    #JM[CartesianIndex.(ip.qN+2:2*ip.qN+1, 1:ip.qN+1)] += ip.nlist/ip.a
    #third set with respect to θ^s (no θ^c derivs)
    #JM[CartesianIndex.(2*ip.qN+2:3*ip.qN+2, 3*ip.qN+3:4*qN+2)] += ip.nlist/ip.a
    #fourth set with respect to θ^c (no θ^s derivs)
    #JM[CartesianIndex.(3*ip.qN+3:4*qN+2, 2*ip.qN+2:3*ip.qN+2)] += ip.nlist/ip.a

    #derivatives w.r.t ν, only δS/δθ (i.e. r^c and r^s) have a non-zero derivative here.
    #not a good name tbh
    #also unsure if there should be some scale factors etc.
    #could almost replace these in the dummy_ν part of drdot, that could be pushing it too far though
    #that would probably add some clarity though@

    #derivative w.r.t r^c
end


"""

Function that packs the CoefficientsT struct into a 1d array for solving.
"""
function new_pack_coeffs!(x::Array{Float64}, coefs::CoefficientsT)

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
    #x .= [coefs.nv ; coefs.rcos ; coefs.θsin[2:end] ; coefs.rsin[2:end] ; coefs.θcos] .+ 1
    #new packing
    x .= [coefs.rcos ; coefs.rsin[2:end] ; coefs.θcos ; coefs.θsin[2:end] ; coefs.nv] .+ 1

end


"""

Function that unpacks the 1d list of coefficients into the CoefficientsT struct.
"""
function new_unpack_coeffs!(x::Array{Float64}, coefs::CoefficientsT, N::Int64)

    #I think we will still want to change the order of this.
    #just need to be a wee bit careful.
    #should be easier with this function and the struct.
    #x[1] = CT.nv[1]
    
    #changed the order hopefully not a mistake!
    
    coefs.rcos[:] = x[1:N] .- 1
    coefs.rsin[2:end] = x[N+1:2*N-1] .- 1
    coefs.θcos[:] = x[2*N-1+1:3*N-1] .-1
    coefs.θsin[2:end] = x[3*N-1+1:4*N-1-1] .-1
    coefs.nv[1] = x[end] - 1
    coefs.rsin[1] = 0.0
    coefs.θsin[1] = 0.0
    
    

    #og packing
    
    #=
    qN = N - 1
    coefs.nv[1] = x[1] - 1
    coefs.rcos .= x[2:qN + 2]  .- 1

    coefs.θsin .= [[0] ; x[qN + 3:2 * qN + 2] .- 1]
    coefs.rsin .= [[0] ; x[2*qN + 3 : 3 * qN + 2] .- 1]
    coefs.θcos .= x[3 * qN + 3 : end] .- 1
    =#
    
end
