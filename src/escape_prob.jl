export k_antibodies_extinction, empirical_log_ext, k_antibodies_extinction_penalized


# since ab_profile is measured in transitions, actually need ne to be 1/2 θ_ts
# θ_ts = 2 n/λ * (β *mu) = 2 n mu /2 = n mu
# fitness: σf = (2 * ne) .* f = θ_ts f  => f = β0/(β mu) =>
# This implies f = should be 2 β_0 /λ mu *
function k_antibodies_extinction(par::PopRates, k, ab_profile; samples = 10^5)
    f = mean_fit(par)
    θ_ref = pop_θ(par)
    λ = par.λ
    rf = (1 - (λ - f)/(λ + f)) * par.κ / θ_ref # rf = σf / θ_ref
    k_antibodies_extinction( θ_ref, k, rf, ab_profile; samples = samples)
end


function k_antibodies_extinction_penalized(par::PopRates, k, ab_profile; samples = 10^5)
    f = mean_fit(par)
    θ_ref = pop_θ(par)
    λ = par.λ
    rf = (1 - (λ - f)/(λ + f)) * par.κ / θ_ref # rf = σf / θ_ref
    Δ_beta = (1/f - 1 /(λ+f)) * (par.ν * λ) / 2 #
    sampler_array = map(x->WrightSampler(x...),
        hcat([θ_ref.* ab_profile for  _ in 1:k]...))
    fit_array = map(x->x[3],
        hcat([θ_ref.* ab_profile for  _ in 1:k]...)) #return an array of fitnesses
    mean( exp(- θ_ref * rf * penalized_escape_freq(sample!.(sampler_array), fit_array ; Δ_beta = Δ_beta)) for ii in 1:samples)
end

function penalized_escape_freq(freq_array, fit_array; Δ_beta = 0) # penalize according to fitness loss
    prod(
        prod(freq_array[ii] .+ (1 .- freq_array[ii]) .* exp.(-Δ_beta .* fit_array[ii]) ) - 
        prod(freq_array[ii]) 
        for ii in eachindex(freq_array))
end

function prob_of_extinction(σf, sampler_array, samples)
    mean( exp(- σf * escape_freq(sample!.(sampler_array))) for ii in 1:samples)
end

#  θ_ref is the base diversity parameter measured in the reference temporal units.
# rf is the ratio of the bare fitness (growth rate) to the mutation rate

function k_antibodies_extinction(θ_ref,k,rf,ab_profile; samples = 10^5)
    sampler_array = map(x->WrightSampler(x...),
        hcat([θ_ref.* ab_profile for  _ in 1:k]...))
    prob_of_extinction(θ_ref .* rf, sampler_array, samples)
end


function escape_freq(freq_array)
    prod(1 .- prod(freq_array, dims = 1))
end


## function of simulations using each individuals exact growth rate

function empirical_log_ext(pop::PopState, par::PopRates;
    sites = [[1,2],[3]] # Example encoding: must escape from and(or(1,2),3)
    )
    escape_fn(genome::BitArray) =  Bool(prod( 
        1 - prod(1-genome[ii] for ii in jj) # 1 if one of the sites is escaped
        for jj in sites)) # 1 if all the sites have atleast 1 escape
    survivors = filter(escape_fn, pop.individuals)
    birth_rates = (par.λ .+ par.f.(survivors)) ./ 2
    death_rates = (par.λ .- par.f.(survivors)) ./ 2
    log_ext = sum(log.(death_rates) .- log.(birth_rates))
end


# Only applicable to a single site, after simulation
# function df_log_ext(df, par::PopRates, fit;
#     sites = 1 # must escape from and(or(1,2),3)
#     )
#     out = log(1 - 2 fit / par.λ) * df.pop_size[] .* df.freq[ind]
# end