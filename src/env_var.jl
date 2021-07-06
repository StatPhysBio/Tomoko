export run_sim, run_sim_total_recombine

function selection_flip(par::PopRates; var_sites = [], prob = 0.1)
    # Generates a new parameter struct with variable sites flipped at random with probability = prob
    n = length(var_sites) # length of vectors with selection
    k = rand(Binomial(n,prob))
    inds = var_sites[StatsBase.sample(1:length(var_sites),k, replace=false)]
    β1 = β1_vector(par.βf)
    β1[inds] .*= -1
    return PopRates(par, β1 = β1)
end

function run_sim(pop::PopState, par::PopRates, timepoints;
	var_sites = [], # vector of indices that are flipping
	flip_prob = 0.0) # probability that they flip at every time point.
    df = initialize_stats() # empty data frame for holding the stats
    for t in timepoints
        run_until!(pop,par,t) # evolve the poulation state until pop.time > t.
        record_stats!(df, pop, par)  # record the statistics in the data_frame
        par = selection_flip(par; var_sites = var_sites, prob = flip_prob) # possibly stochastically change the environment.
    end
    return df # return the df, now heavy with juicy statistics for you to plot or analyze
end

# default definition when pop is provided
function run_sim(par::PopRates, timepoints; kwdargs...) # population loci set to wright equilibrium, 
	pop = initialize_pop(par)
	run_sim(pop, par, timepoints; kwdargs...)
end


"""
bottleneck removes all individuals for which survives = false
"""
function bottleneck!(pop::PopState, survives::Function)
    pop.individuals = pop.individuals[survives.(pop.individuals)]
end

function shuffle_col!(data)
    vdata = Align(data,1)
    for i = 1:size(vdata, 1)
        shuffle!(@view vdata[i,:])
    end
    return
end

function total_recombine(pop::PopState)
    shuffle_col!(pop.individuals)
end

function run_sim_total_recombine(par::PopRates, timepoints; kwdargs...) # population loci set to wright equilibrium, 
	pop = initialize_pop(par)
	run_sim_total_recombine(pop, par, timepoints; kwdargs...)
end

function run_sim_total_recombine(pop::PopState, par::PopRates, timepoints;
	var_sites = [], # vector of indices that are flipping
	flip_prob = 0.0) # probability that they flip at every time point.
    df = initialize_stats() # empty data frame for holding the stats
    for t in timepoints
        run_until!(pop,par,t) # evolve the poulation state until pop.time > t.
        record_stats!(df, pop, par)  # record the statistics in the data_frame
        total_recombine(pop) # get near limit of infinite recombination.
        par = selection_flip(par; var_sites = var_sites, prob = flip_prob) # possibly stochastically change the environment.
    end
    return df # return the df, now heavy with juicy statistics for you to plot or analyze
end