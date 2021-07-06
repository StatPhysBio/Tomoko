export estimate_rs

"""
This is a MLE for fitting the selection parameters to a time trace, given the neutral parameters
"""

##
"""
The log-likelihood difference is computed using multiple importance sampling 
with a mixture of the neutral-beta and the binomial-reweighted neutral sample
1/2 * B(Dq, D(1-q)) + 1/2 * B(Dq+k, D(1-q) + N-k).
This give the benefit of using the same samples for both likelihoods, greatly reducing the noise.

"""

logbinomial(j, k) =  loggamma(1+j+k) - loggamma(1+j) - loggamma(1+k)
##
function randlogGamma(a,n;scale = 1)
    if a > .2
        return log.(rand(Gamma(a),n))
    else
        L = 1 / a - 1
        w = a / MathConstants.e / (1 - a)
        ww = 1 / (1 + w)
        η = z -> ( z >= 0 ? exp(-z) : w * L * exp(L * z))
        h = z -> exp(-z - exp(-z / a))
        function rh(a)
            while true
                U = rand()
                z = (U <= ww) ? -log(U / ww) : log(rand()) / L
                h(z) / η(z) > rand() && return(z)
            end
        end
        Z = Array{Float64}(undef,n)
        for ii in eachindex(Z)
            Z[ii] = log(scale) - rh(a) / a
        end 
        return Z
    end
end

# function sigmoid(x) # numerically stable sigmoid
#     if x >= 0
#         1 / ( 1 + exp(-x) )
#     else
#         exp(x) / ( 1 + exp(x) )
#     end
# end

# Likelihood parameterized strictly in terms of exponents
function σ_likelihood(θ_wt, θ_mut, k_wt, k_mut) # mutant and wildtype counts
    n_samp = 1000 # hyperparameter 
    loggam_a = vcat(
        randlogGamma(θ_wt, n_samp),
        randlogGamma(θ_wt+k_wt, n_samp)
    )
    loggam_b = vcat(
        randlogGamma(θ_mut, n_samp),
        randlogGamma(θ_mut+k_mut, n_samp)
    )
    xs = 1 ./(1 .+ exp.(loggam_b .- loggam_a))
    δ_z = logbeta(θ_wt, θ_mut) - logbeta(θ_wt+k_wt, θ_mut +k_mut)
    δ_ls = map((a,b) -> k_wt*a + k_mut*b - (k_mut + k_wt)*(max(a,b) + log1p(exp(-abs(a-b)))) + δ_z, loggam_a, loggam_b) # log(p(d+k)/p(d)) 
    σ -> log(sum( exp(x*σ) / (1 + exp(δ_l)) for (x,δ_l) in zip(xs,δ_ls))) -
        log(sum( exp(x*σ) / (1 + exp(-δ_l)) for (x,δ_l) in zip(xs,δ_ls))) - 
        logbinomial(k_wt,k_mut)
end



# defines an importance sampler structure for defining a likelihood

struct LikelihoodSample
    w_θ::Array{Float64,1}
    w_k::Array{Float64,1}
    x::Array{Float64,1}
    θ_ts::Float64
end

function LikelihoodSample(μ_wt, μ_mut, θ_ts, k_wt, k_mut; n_samp = 1000)
    (θ_wt, θ_mut) = θ_ts .* (μ_wt, μ_mut)
    # quenched and variance reduced sampler
    loggam_a = vcat(
        randlogGamma(θ_wt, n_samp),
        randlogGamma(θ_wt+k_wt, n_samp)
    )
    loggam_b = vcat(
        randlogGamma(θ_mut, n_samp),
        randlogGamma(θ_mut+k_mut, n_samp)
    )
    δ_z = logbeta(θ_wt, θ_mut) - logbeta(θ_wt+k_wt, θ_mut +k_mut)
    δs = map((a,b) -> k_wt*a + k_mut*b - (k_mut + k_wt)*(max(a,b) + log1p(exp(-abs(a-b)))) + δ_z, 
        loggam_a, loggam_b) # log(p(d+k)/p(d)) 
    LikelihoodSample( 
        2 .* sigmoid.(-δs) , 
        2 .* sigmoid.(δs) ,
        sigmoid.(loggam_a .- loggam_b),θ_ts)
end

function like_rs(rs, ls::LikelihoodSample)
    θ_ts = ls.θ_ts
    z_θ = mean(w_θ*exp(θ_ts*rs*x) for (w_θ, x) in zip(ls.w_θ, ls.x))
    z_k = mean(w_k*exp(θ_ts*rs*x) for (w_k, x) in zip(ls.w_k, ls.x))
    return log(z_k/z_θ)
end

function grad_hess_rs(rs, ls::LikelihoodSample)
    θ_ts = ls.θ_ts
    z_θ = mean(w_θ*exp(θ_ts*rs*x) for (w_θ, x) in zip(ls.w_θ, ls.x))
    z_k = mean(w_k*exp(θ_ts*rs*x) for (w_k, x) in zip(ls.w_k, ls.x))
    x_θ = mean(w_θ*exp(θ_ts*rs*x)*x for (w_θ, x) in zip(ls.w_θ, ls.x)) / z_θ
    x_k = mean(w_k*exp(θ_ts*rs*x)*x for (w_k, x) in zip(ls.w_k, ls.x)) / z_k
    xx_θ = mean(w_θ*exp(θ_ts*rs*x)*x^2 for (w_θ, x) in zip(ls.w_θ, ls.x)) / z_θ
    xx_k = mean(w_k*exp(θ_ts*rs*x)*x^2 for (w_k, x) in zip(ls.w_k, ls.x)) / z_k
    g = (x_k - x_θ)*θ_ts
    h = ((xx_k - x_k^2) - (xx_θ - x_θ^2)) *θ_ts^2
    return (g,h)
end

function grad_hess_pop(rs,ls::Array{Array{LikelihoodSample,1},1}; average = true)
    g = 0
    h = 0
    for lsp in ls
        gp = 0
        hp = 0
        for lspt in lsp
            (gpt,hpt) =  grad_hess_rs(rs,lspt)
            gp += gpt
            hp += hpt
        end
        if average
            g += gp / length(lsp) #
            h += hp /length(lsp) # time averaging can be turned off here by commenting "/length(lsp)"
        else
            g += gp  #
            h += hp  # time averaging can be turned off here by commenting "/length(lsp)"
        end
    end
    return (g,h)
end


function rs_newton_search(ls; average = true)
    rs = 0
    while true
        (g,h) = grad_hess_pop(rs,ls ; average = average)
        if |(g^2/abs(h) < 10^(-8), abs(rs)>10^4)
            return (rs, g, h)
        end
        rs -= g/h
    end
end

function like_rs_pop(rs, ls; avg = 0) 
    # ls is an array of arrays of likelihood sampers, [patients[timpoints]], where each time point 
    # avg is a number between 0 and 1
    l = 0
    for lsp in ls
        lp = 0
        for lspt in lsp
            lp += like_rs(rs,lspt)
        end
        l += ( 1/length(lsp) * avg + 1 - avg)*lp
    end
    return l
end

function baysian_pop_rs(ls; samples = 10^4, burn_in = 10^3, avg = 0) # ls is an array of likelihood sampers
    ll = rs -> like_rs_pop(rs, ls; avg = avg)
    baysian_sample_rs(ll; samples = samples, burn_in = burn_in)
end

function baysian_sample_rs(ll; samples = 10^4, burn_in = 50)
    # takes a log-likelihood and returns the sampled distribuiton with
    run  = Float64[]
    rs = 100
    for ii in 1:(samples + burn_in)
        rs = mh_step(ll,rs, σ = 50)
        push!(run,rs)
    end
    return run[burn_in+1:end]
end

function mh_step(ll, θ; σ = 50) # 1-d mh_stepper
    oldll = ll(θ)
    new = θ + randn().* σ
    newll = ll(new)
    rand() < exp(newll - oldll) ? out = new : out = θ # ll difference, always accept if positive
    return out
end


# function get_countdata(df, locus, time_points)
#     ind = map( t->findfirst(x -> x>t, df.time), time_points)
#     ind = ind[ind .!= nothing]
#     Ds = df.D[ind] # population diversity
#     # subsampling element (this may be neccesary to limit the large number of genomes)
#     Ns = df.pop_size[ind] # total number of indivudals
#     Ks = round.(Int,map(x -> x[locus], df.freq[ind]) .* Ns) #number of mutants.
#     Js = Ns .- Ks
#     return (Ds, Js, Ks)
# end
    


function make_likelihood(df; locus = 1, timepoints = 100:100:1000, sub_n = 100)
    # Return the likelihood function for a single time trace
    ind = map( t->findfirst(x -> x>t, df.time),timepoints)
    ind = ind[ind .!= nothing]
    if length(ind) <= 1
        return (x->0)
    end
    θs = df.θ[ind] # population diversity
    # subsampling element (respresenting incomplete information in a typical genetic data point)
    Ns = df.pop_size[ind] # total number of indivudals
    Ks = round.(Int,map(x -> x[locus], df.freq[ind]) .* Ns) #number of mutants.
    Js = Ns .- Ks
    ks = map((k,j) -> rand(Hypergeometric(k, j, min(sub_n,k+j))) .* [1,-1] .+ [0,  min(sub_n,k+j)] , Ks, Js)
    # if pop.size is smaller than the subsample size (sub_n), than default to pop.size
    fn_list = map((θ_ref,k)->rs_likelihood(θ_ref, [1,1], k...), θs, ks)
    fn_normalize = rs_likelihood(10^-5, [1,1], 1,1)
    # this function creates a minimum at finite Δ even if there is no mutants observed
    # D sets rough scale at which fitness effect will be invisible.
    rs -> sum(f(rs) for f in fn_list) + 0.01*fn_normalize(rs)
end

function estimate_rs(dfs::Vector; locus = 1, timepoints = 100:100:1000, sub_n = 100)
    fnlist = map(x-> make_likelihood(x; locus = locus, timepoints = timepoints, sub_n = sub_n),dfs)
    like(rs) = mapreduce(f->f(rs...), +, fnlist)
    s0 = [0.0]
    y = Optim.minimizer(optimize(like, s0, LBFGS(), autodiff = :forward))
    return -y[1]
end

function get_Δ(dfs::Vector; locus = 1, timepoints = 100:100:1000, sub_n = 100)
    fnlist = map(x-> make_likelihood(x; locus = locus, timepoints = timepoints, sub_n = sub_n),dfs)
    like(Δ) = mapreduce(f->f(Δ...), +, fnlist)
    return like
end