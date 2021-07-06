export pop_stats # global variable listing the a
export record_stats!, initialize_stats

"""
To extend the statistics to include a function, the function must live inside the Tomoko module.
This can be done either by modifying this file or in the REPL using the following interface:

```
julia> Tomoko.eval(:(
function est_fitnesses(pop::PopState, par::PopRates)
    est_fit(pop::PopState, par::PopRates, sites = sites_of_interest)
end
))

julia> push!(pop_stats,:new_stat)
```

at which point record_stats! will start recording another column in new data-frames
"""

## Statistics on populations
# used to construct a dataframe of simulation results.
global pop_stats = [:time, :pop_size, :freq, :θ, :mean_fit, :var_fit, :mut_rate, :noise]
# make this available to extend with one's own functions.

initialize_stats() = DataFrame();

function record_stats!(df::DataFrame, pop::PopState, par::PopRates,stats = pop_stats)
    push!(df,
    (;[(f,getfield(Tomoko,f)(pop,par)) for f in stats]...) # construct a named tuple
    )
end

function time(pop::PopState, par::PopRates)
    pop.time
end

function pop_size(pop::PopState, par::PopRates)
    pop.size
end

"""
If fitnesses are known, we calculate diversity only from the neutral sites.
"""

function θ(pop::PopState, par::PopRates; 
    k = length(pop.individuals)) # size of population to take counts one, defaults to all
    stat_D(pop, sites = par.βf.nind, k = k)/2
end

function mut_rate(pop::PopState, par::PopRates)
    par.ν * par.λ / 2
end

function noise(pop::PopState, par::PopRates)
    par.λ
end

function freq(pop::PopState, par::PopRates)
    sum(pop.individuals)/pop.size
end

function mean_fit(pop::PopState, par::PopRates)
    mean(pop.fitness)
end

function var_fit(pop::PopState,par::PopRates)
    var(pop.fitness)
end


"""
The MLE fit for the beta-binomial distribution

This follows T. Minka's fixed point algorithm for fitting Polya distributions.
"""

function neutral_update_D(data::Array{T,2}, s::Real) where T
    # data has size = (number of categories, number of samples)
    mk = [0.5; 0.5] 
    ns = map(sum, Slices(data, 1))
    f1 = sum(digamma(s) - digamma(ns[ii] + s) +
            sum(mk[jj] * ( digamma(data[jj,ii] + s*mk[jj]) - digamma(s*mk[jj]))
            for jj in eachindex(mk))  
        for ii in eachindex(ns))
    f2 = sum(polygamma(1,s) - polygamma(1,ns[ii] + s) +
            sum(mk[jj]^2 * ( polygamma(1, data[jj,ii] + s*mk[jj]) - polygamma(1,s*mk[jj]))
                for jj in eachindex(mk)) for ii in eachindex(ns))
    a = -s^2*f2
    c = f1 - a/s
    if c >= 0
        a = s^3*(s*f2 + 2*f1)
        c = -(s^2 *f1 + a/s)
    end
    return -a/c
end

function neutral_D(data::Array{T,2}) where T
    s = 1
    while true
        s_new = neutral_update_D(data, s)
        abs(s - s_new) < 10^-9 ? break : s = s_new
    end
    return s
end

function stat_D(pop::PopState; 
    sites = 1:length(pop.individuals[1]), # vector of sites to use in likelihood 
    k = length(pop.individuals)) # size of population to take counts one, defaults to all
        vec = sum(StatsBase.sample(pop.individuals, k, replace=false))[sites]
        vec = permutedims(hcat(vec , k .- vec))
        neutral_D(vec)
end


##
function estimate_θ(n_sites,n_genomes,true_θ)
    results = Float64[]
    for ii in 1:1000
    data = mapreduce(hcat,rand(Beta(true_θ,true_θ),n_sites)) do p
        k = rand(Binomial(n_genomes,p))
        [k,n_genomes-k]
    end
    push!(results,Tomoko.neutral_D(data))
    end
    return results
end