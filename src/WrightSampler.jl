#=
This document implements a WrightSampler mutable struct type.
This object is used to generate Gibbs-sampler output on Wright equilibria
Also included are methods of the current state to generate stochastic gradients
And fisher information estimates
 =#


# This lists the neccesary parameters 
# for sampling from the Wright-equilibrium distribution. The partition function is the most costly to compute and so it cahched
# θ_mut, θ_wt, σ are forward and backward diversity parameters respectively. At the moment, there is only code for  σ>0.
mutable struct WrightSampler
    θ_wt::Float64
    θ_mut::Float64
    σ::Float64
    k::Int64   # auxiliary variable for gibbs sampling
    x::Float64 # The freqeuncy of wt
end


function WrightSampler(θ_wt, θ_mut, σ)
#    z = gamma(θ_wt) * HypergeometricFunctions.pFq([θ_wt],[θ_mut + θ_wt],σ) /gamma(θ_wt +θ_mut)
    ws = WrightSampler(θ_wt, θ_mut, σ, 0, .5)
    for ii in 1:10
    incrementk!(ws)
    incrementx!(ws)
    end
    return ws
end

# This sampler is robuset for small-exponent gamma distributions
# Code was copied from somewhere unknown (citation needed)
function randlogGamma(a;scale = 1)
    if a > .2
        return log.(rand(Gamma(a)))
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
        return log(scale) - rh(a) / a
    end
end

function sigmoid(x) # numerically stable sigmoid
    if x >= 0
        1 / ( 1 + exp(-x) )
    else
        exp(x) / ( 1 + exp(x) )
    end
end

function incrementk!(ws::WrightSampler)
    if ws.σ >= 0 
        ws.k = rand(Poisson(abs(ws.σ * ws.x)))
    else
        ws.k = rand(Poisson(abs(ws.σ * (1-ws.x))))
    end
end

function incrementx!(ws::WrightSampler)
    if ws.σ >= 0 
        a = randlogGamma(ws.θ_wt + ws.k;scale = 1)
        b = randlogGamma(ws.θ_mut;scale = 1)
    else
        a = randlogGamma(ws.θ_wt;scale = 1)
        b = randlogGamma(ws.θ_mut + ws.k;scale = 1)
    end
    ws.x = sigmoid(a-b)
end

function sample!(ws::WrightSampler)
    incrementk!(ws)
    incrementx!(ws)
    return ws.x
end


function rejection_sample_k(θ_wt, θ_mut,σ)
    k = 0 
    while true
        k = rand(Poisson(σ))
        k == 0 ? break : nothing
        p = exp( loggamma(θ_wt + θ_mut) + loggamma(k + θ_wt) - loggamma(k + θ_wt + θ_mut) -loggamma(θ_wt))
        rand(Bernoulli(p)) ? break : nothing
    end
    return k
end

## Generate rebound values

struct MutantState
	n::Int # number of individuals
	t::Float64
end

function reboundtime_step_closure(;theta = 0.02, mut_per_gen = 10^5, extra_mut_flux = 1.0, delta_per_mut = 100, pop_size::Int = 10^4, decay_rate = .5, x0 = 10^-3)
	# x0 is the inital mutant fraction
	let f = 1/3, 
	µ = mut_per_gen * f, 
	θ = .02, K = pop_size, 
	λ = μ*2*K/θ, r =decay_rate, 
	delta = delta_per_mut*μ # paramters in absolute units (running on f = 1/3 days)
	if λ<f 
		error("pop_size too small to support growth rate") 
	end
	function step_state(state::MutantState)
		mut = extra_mut_flux*K*μ*exp(-r*state.t)
		rate = mut + λ*state.n
		delta_t = randexp()/rate
		if mut/rate > rand()
			return MutantState(state.n+1,state.t+delta_t)
		else
			birth = (λ + (f - delta) - 
					((1-x0)*exp(-r*state.t)*f + state.n*(f - delta)/K)
					)/2
			if birth/λ > rand()
				return MutantState(state.n+1,state.t+delta_t)
			else
				return MutantState(state.n-1,state.t+delta_t)
			end
		end
	end
		
	end
	
end

function run_reboundtimes(start::Int ; theta = 0.02, nsample = 1000, mut_per_gen = 10^(-5.), extra_mut_flux = 1.0, delta_per_mut = 100, pop_size=10^4, decay_rate = .5)
	steppy = reboundtime_step_closure(;theta, mut_per_gen, extra_mut_flux, delta_per_mut, pop_size, decay_rate, x0 = start/pop_size)
    results = Float64[]
    q = 0.5
	for _ in 1:nsample
	state = MutantState(start,0);
	while (&)(state.t < 56, state.n <  pop_size*q)
		state = steppy(state)
	end
	if state.t < 56
		t_reb = 3*log(1 + exp(state.t / 3)*(1-q)/q) # Interpolate to deterministic rebound time
	else
		t_reb = 57
	end
	push!(results,t_reb)
	end
	return results
end



# function ab_escapes(theta::Array{Float64,1}, abprofiles::Array{Array{Float64,1},1}...; samples = 1)
#     out = Float64[]
#     for θ in theta
#         samplingarray = [map(x->WrightSampler( (θ.*x)... ), ab) for ab in abprofiles] # array of array of wright samplers
#         for _ in 1:samples
#         push!(out, prod( sample!(ab)) for ab in abprofiles)) # escape probability is the probability of being escaped to all of the abs
#         #each ab consists of an array of arrays of site profiles 
#         end
#     end
# end
