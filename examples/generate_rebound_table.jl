using Tomoko
using HDF5
using Statistics
using BenchmarkTools
using BitBasis
using Tomoko
using StatsBase
using Random

## set of helper functions and path definitions for reading in data



snpanalysis = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/snpanalysis.h5"
trialpath(x) = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialanalysis$x.h5"
trial_list = ["10-1074","3BNC","combo"]

h5open(snpanalysis, "r") do fid
	global ablist = keys(fid)
	ablist = [ab for ab in ablist if (&)(ab != "101074", ab != "VRC01_dms")]
end 


function get_start_theta(trial)
	if trial != "all"
    h5trial = h5open(trialpath(trial), "r")
    thetas = Float64[]
    for patient in h5trial["genetic"]
        for day in patient
            if (read(day["day"]) < 2) & haskey(day, "theta")
                push!(thetas,read(day["theta"]))
            end
        end
    end
    close(h5trial)
    return thetas
	else 
	return mapreduce(get_start_theta, vcat, trial_list)
	end
end

function get_antibody_profile_mle(ab; selection_multiplier = 1.0)
    # return the mle parameters for the antibody
    sites = Array{Float64,1}[]
    h5open(snpanalysis, "r") do fid
            for site in fid[ab]
                rs = read(site["rsel"])[1] * selection_multiplier
                mut = reverse(read(site["mut"])) # [backmut, forwardmut]
                if rs > 0
                    push!(sites,vcat(mut...,rs))
                end
            end
    end
    return sites
end

function get_antibody_profile_bayes(ab; selection_multiplier = 1.0)
    # return the mle parameters for the antibody
    sites = Array{Float64,1}[]
    h5open(snpanalysis, "r") do fid
            for site in fid[ab]
                rs = rand(read(site["bayes_avg1.0"])) * selection_multiplier
                mut = reverse(read(site["mut"])) # [backmut, forwardmut]
                if rs > 0
                    push!(sites,vcat(mut...,rs))
                end
            end
    end
    return sites
end


function get_observed_rebound_times(trial)
    # return the mle parameters for the antibody
	t_rebounds = Float64[]
    h5trial = h5open(trialpath(trial), "r")
	growth_rate = read(h5trial["viraemic/growth"])
    for patient_key in keys(h5trial["viraemic/patients"])
		patient = h5trial["viraemic/patients/$patient_key"]
		if (|)(trial != "3BNC", !(in(["2A1","2A3","2A4"])(patient_key))) 
			push!(t_rebounds,-log(read(patient["x"])) / growth_rate)
		end
    end
    close(h5trial)
    return t_rebounds
end

function wright_sampler_array(θ, ab_list; selection_multiplier = 1)
	out = Tomoko.WrightSampler[]
    for ab in ab_list 
		for site in get_antibody_profile_mle(ab; selection_multiplier)
			push!(out,Tomoko.WrightSampler((θ.* site)...))
		end
	end
	return out
end

function wright_sampler_array_bayes_closure(ab_list; selection_multiplier = 1)
	sites = [] # flat vector of samplers for each site
    for ab in ab_list 
		for site in get_antibody_profile_bayes(ab; selection_multiplier)
			push!(sites,site)
		end
	end
	return θ -> map(site -> Tomoko.WrightSampler((θ.* site)...), sites)
end


##
#-------------------------------
# Actual gillespie code 
#-------------------------------


struct Virus
	genotype::Int64 # this is the genotype of bits.
	fitness::Float64 # the additive fitness of the virus
	escape::Bool # whether or not the virus is escaped
end

mutable struct ViralPopulation
	time::Float64 # current time
    fitness_total::Float64 # current fitness. Divided by cap gives φ
    pop::Array{Virus,1} # State of population length is the total population size
    capacity::Int 	# carrying capacity
    r::Float64 	# rate of decay of neutralized viruses
    λ::Float64 	# noise of system
	fitness_function::Function
	mutation_operator::Function
	escape_function::Function
end


function mutation_operator_closure(ab_list, ν)
    let mut::Array{Float64,1} =  vcat([vcat(map(x->x[1:2], 
					get_antibody_profile_mle(ab))...) for ab in ab_list]...),
		l = Int64(length(mut)/2), tot = sum(mut)
		# this is a let-block to capture a fixed
		function mutate!(vp::ViralPopulation) # attempt to mutate a virus at index ind
            if rand() < tot*ν # if there is a shot at mutation
                ind = rand(1:length(vp.pop)) # choose a random virus
				v = vp.pop[ind]
				mutvec = mut[ (1:2:2*l) .+ bitarray(v.genotype,l)] # get the mutational parameters
				if rand() < sum(mutvec)/tot # if a mutation actually happens
					# pick a genotype to invert 
					ii = sample(Weights(mutvec))
					# flippy the dippy and make a new virus
					new_genotype = flip(v.genotype, 2^(ii-1))
					new_fitness = vp.fitness_function(new_genotype)
                    new_escape = vp.escape_function(new_genotype)
                    vp.fitness_total += new_fitness - v.fitness
                    vp.pop[ind] = Virus(new_genotype,new_fitness,new_escape)
				end
			end
        end
        return mutate!
	end
end

function escape_function_closure(ab_list)
	# true or false on being escaped or not
	# determines binding in simulations
	ab_site_pattern = [length(get_antibody_profile_mle(ab)) for ab in ab_list]
	let bitpatterns = Int64[packbits(vcat([
							if (n == m) 
								ones(ab_site_pattern[m]) else
								zeros(ab_site_pattern[m]) 
							end
		for m in 1:length(ab_site_pattern)]...))
				for n in 1:length(ab_site_pattern)]
		# construct the bit masks for each antibody
		genotype -> mapreduce( x->anyone(genotype,x), (&), bitpatterns)
	end
end

function fitness_function_closure(ab_list, f, μ; selection_multiplier = 1.0)
	# absolute fitness and mutation rate
    let deltavec::Array{Float64,1} =  vcat([map(x->x[3], get_antibody_profile_mle(ab; selection_multiplier)) for ab in ab_list]...) .* μ
		# this is a let-block to capture a fixed parameter
		function fitness(b::Int64)
			result = 0.0
				for ii in baddrs(b)
				result += deltavec[ii]
				end
            return f - result
        end
        return fitness
	end
end

function standing_variation(θ, ab; selection_multiplier = 1)
    map(x-> random_wright_sample((θ .* x )...), get_antibody_profile_mle(ab; selection_multiplier))
end

function random_genotype(frequency_vector::Vector{Float64})
	x::Int64 = 0
	for ii in Iterators.reverse(frequency_vector) #start from the most significant bit
		x += rand() > ii
		x <<= 1 # bitshift to a larger power of 2
	end
	x >>=1 #undo the last bitshift
	return x
end

function initialize_viral_population(θ::Float64, ab_list;
		mut_per_gen = 3*10^(-6), # the mutation rate measured in growth rate 10^-6
		decayrate = 0.4,
		f = 1.0/3, selection_multiplier = 1.0)
	μ = mut_per_gen * f # absolute transition mutation rate in equilirbrium
	λ = 2 * decayrate  #  λ = μ*2*capacity/θ   =>    λ * θ/ (2 μ) = capacity
	capacity = ceil(Int64,λ * θ / (2 * μ))
    ν = θ / capacity / 2# per-birth-event mutation rate; ν = μ/λ = θ / 2 capacity
	if λ<f 
		error("pop_size too small to support growth rate") 
	end
	
	start_freq = vcat([standing_variation(θ, ab ; selection_multiplier) for ab in ab_list]...)
	fitness_function = fitness_function_closure(ab_list, f, μ ; selection_multiplier)
	mutation_operator = mutation_operator_closure(ab_list, ν)
	escape_function = escape_function_closure(ab_list)
	pop = Virus[]
	genotype_length = length(start_freq) 
	basis = 2 .^ collect(0:(genotype_length-1))
	for ii in 1:capacity
		genotype = random_genotype(start_freq)
		fitness = fitness_function(genotype)
		escape = escape_function(genotype)
		push!(pop, Virus(genotype, fitness, escape))
	end
	fitness_total = sum(virus.fitness for virus in pop)
	ViralPopulation(0, fitness_total, pop, capacity, decayrate, λ, fitness_function, mutation_operator, escape_function)
end

function restart_viral_population(original_vp::ViralPopulation, start_freq::Vector{Float64})
	# this looks like it is still based on the static struct form
	# really only need to update time and pop state
	# have to do some profiling so see if this is really a problem
	# fair amount of memory allocation,
fitness_function = original_vp.fitness_function
mutation_operator = original_vp.mutation_operator
escape_function = original_vp.escape_function
capacity = original_vp.capacity
λ = original_vp.λ
decayrate = original_vp.r
pop = Vector{Virus}(undef,capacity)
genotype_length = length(start_freq) 
basis = 2 .^ collect(0:(genotype_length-1))
for ii in 1:capacity
	genotype = random_genotype(start_freq)
	fitness = fitness_function(genotype)
	escape = escape_function(genotype)
	pop[ii] = Virus(genotype, fitness, escape)
end
fitness_total = sum(virus.fitness for virus in pop)
ViralPopulation(0, fitness_total, pop, capacity, decayrate, λ, fitness_function, mutation_operator, escape_function)
end

function random_swap!(array::Vector) # bring an element of a vector to the end
    # My thinking is that this will be good for memory management
    ii = rand(1:length(array))
    @inbounds (array[end],array[ii]) = (array[ii],array[end]);
end

function gillespie_step!(vp::ViralPopulation)
    N = length(vp.pop)
    vp.time += randexp() / (N*vp.λ) # step the population forward a little bit
if N > 0
	vp.mutation_operator(vp) # see if you can mutate something
    random_swap!(vp.pop)
    virus = pop!(vp.pop) # get the virus stack-style
	if virus.escape
    	if rand() < ((vp.λ - vp.fitness_total/vp.capacity + virus.fitness) / vp.λ / 2 )
        	push!(vp.pop,virus)
        	push!(vp.pop,virus)
			vp.fitness_total += virus.fitness
		else
			vp.fitness_total -= virus.fitness
    	end # else get rid of that shit, it is dead
	else # if the virus is neutralized
		if rand() >  vp.r /vp.λ # check and see if you are actually getting rid of it at the proper rate
			push!(vp.pop,virus) # put that virus back if you tried to kill it too quickly
		else
			vp.fitness_total -= virus.fitness
    	end # else get rid of that shit, it is dead
	end
end
end

function run_viruses_stop(vp::ViralPopulation, end_time)
	while vp.time < end_time
		gillespie_step!(vp)
	end
end

##
# --------------------------------------------------------------
# Funcitons for generating traces and measuring quantities
# --------------------------------------------------------------




function wt_esc_counts(vp)
	esc = 0
	for virus in vp.pop
		if virus.escape
			esc += 1
		end
	end
	[esc, length(vp.pop) - esc]
end

function virus_time_trace(vp::ViralPopulation, time_points; breakpoint = 5.0)
	# breakpoint sets a threshold at which to end the simulation and save on time
	# because escape has already been acieved
	out = []
	for time in time_points
		run_viruses_stop(vp, time)
		(esc,wt) = wt_esc_counts(vp)
		push!(out, [esc,wt])
		if esc/vp.capacity > breakpoint
			break
		end
	end
	return hcat(out...)
end

function simulate_viral_rebound_time(vp::ViralPopulation; breakpoint = 0.8)
	# breakpoint is a keyword that determines when a simulation should be treated deterministically
	# it is set very conservatively here. The main purpose is to prevent running simulations 
	# at full capacity all the way to the endpoint time (56 days)
	(esc,wt) = wt_esc_counts(vp)
	trace = virus_time_trace(vp, 0:56; breakpoint)
	if esc/vp.capacity > breakpoint
		t_rebound = 0.0
	elseif trace[1,end]/vp.capacity < breakpoint
		t_rebound = 57
	else
		ind = findfirst(col->col[1]/vp.capacity > breakpoint, collect(eachcol(trace)))
		t_cross = (0:56)[ind]
		q_cross = trace[1,ind]/vp.capacity
		t_rebound = 3*log(1 + exp(t_cross / 3) * (1-q_cross)/q_cross)
	end 
	return (t_rebound,trace)
end

# function simulate_trial_rebounds(trial, ab_list; samples_per_patient = 10)
# 	theta_vec = get_start_theta(trial)
# 	rebound_times = Float64[]
# 	for θ in theta_vec
# 		for _ in  1:samples_per_patient
# 			vp = initialize_viral_population(θ,ab_list)
# 			(t_rebound,_) =  simulate_viral_rebound_time(vp)
# 			push!(rebound_times,t_rebound)
# 		end
# 	end
# 	return rebound_times
# end

function trial_traces(trial, ab_list; samples_per_patient = 10, 
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0, 
	no_mutations = false, 
	breakpoint = 2,
	timepoints = 0:.5:56
	)
#= returns a 3 d array of viral rebound traces [:,1,:] being wt and [:,2,:] being mutant  =#
	theta_vec = get_start_theta(trial)
	rebound_times = ones(length(theta_vec), samples_per_patient )
	trace = []
	out_mat = permutedims(hcat(collect(timepoints),collect(timepoints)))
	for (ii,θ_0) in enumerate(theta_vec)
		θ = θ_0 * diversity_multiplier
		ws = wright_sampler_array(θ, ab_list; selection_multiplier)
		xxvec = [Tomoko.sample!.(ws) for jj in 1:samples_per_patient ]
		vp = initialize_viral_population(θ,ab_list; selection_multiplier)
		if no_mutations
			vp.mutation_operator = x -> nothing # if mutations are turned off set mutations to do nothing
		end

		for jj in  1:samples_per_patient
			xx = xxvec[jj]
			new = restart_viral_population(vp, xx)
			trace =  virus_time_trace(new, timepoints; breakpoint = 5.0) ./ new.capacity
			out_mat = cat(out_mat, trace, dims = 3)
		end
	end
	return out_mat
end



function simulate_trial_rebounds_matrix(trial, ab_list; samples_per_patient = 10, 
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0, 
	no_mutations = false, 
	breakpoint = .8
	)
#= main function returns a matrix of samples organized by patient
	This allows us to simulate trial outcomes when choosing without replacement in fit analysis.
	But also gets flattened for raw histogram construction.
	no_mutations keyword turns off mutations after treatment.
	multiplier keywords are used in the sensitivity analysis. =#
	theta_vec = get_start_theta(trial)
	rebound_times = ones(length(theta_vec), samples_per_patient )
	for (ii,θ_0) in enumerate(theta_vec)
		θ = θ_0 * diversity_multiplier
		ws = wright_sampler_array(θ, ab_list; selection_multiplier)
		xxvec = [Tomoko.sample!.(ws) for jj in 1:samples_per_patient ]
		vp = initialize_viral_population(θ,ab_list; selection_multiplier)
		if no_mutations
			vp.mutation_operator = x -> nothing # if mutations are turned off set mutations to do nothing
		end
		Threads.@threads for jj in  1:samples_per_patient
			xx = xxvec[jj]
			new = restart_viral_population(vp, xx)
			(t_rebound,_) =  simulate_viral_rebound_time(new; breakpoint)
			rebound_times[ii,jj] = t_rebound
		end
	end
	return rebound_times
end

function simulate_trial_rebounds(trial, ab_list; samples_per_patient = 10, 
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0
	)
	rebound_times = vec(simulate_trial_rebounds_matrix(trial, ab_list;
		samples_per_patient, selection_multiplier, diversity_multiplier))
	return rebound_times
end

function simulate_trial_rebound_prob(trial, ab_list; samples_per_patient = 10, 
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0, 
	no_mutations = false, breakpoint = .05
	)
	rebound_times = vec(simulate_trial_rebounds_matrix(trial, ab_list;
		samples_per_patient, no_mutations, breakpoint))
	return mean(x < 56.0 for x in rebound_times)
end

function simulate_trial_rebound_prob_bayes(ab_list; samples_per_patient = 10, 
	trial = "all",
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0, 
	no_mutations = false, breakpoint = 0.01
	)
	rebound_times = vec(simulate_trial_rebounds_bayes(trial, ab_list; diversity_multiplier,
		samples_per_patient, no_mutations, breakpoint))
	return mean(x < 56.0 for x in rebound_times)
end



function simulate_trial_rebounds_bayes(trial, ab_list; samples_per_patient = 10, 
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0, 
	no_mutations = false, 
	breakpoint = .8
	)
#= main function returns a matrix of samples organized by patient
	This allows us to simulate trial outcomes when choosing without replacement in fit analysis.
	But also gets flattened for raw histogram construction.
	no_mutations keyword turns off mutations after treatment.
	multiplier keywords are used in the sensitivity analysis. =#
	theta_vec = get_start_theta(trial)
	rebound_times = ones(length(theta_vec), samples_per_patient )
	ws_fun = wright_sampler_array_bayes_closure(ab_list; selection_multiplier)
	# creates a random function of θ to generate correlated frequencies
	for (ii,θ_0) in enumerate(theta_vec)
		θ = θ_0 * diversity_multiplier
		ws = ws_fun(θ)
		xxvec = [Tomoko.sample!.(ws) for jj in 1:samples_per_patient ]
		vp = initialize_viral_population(θ,ab_list; selection_multiplier)
		if no_mutations
			vp.mutation_operator = x -> nothing # if mutations are turned off set mutations to do nothing
		end
		Threads.@threads for jj in  1:samples_per_patient
			xx = xxvec[jj]
			new = restart_viral_population(vp, xx)
			(t_rebound,_) =  simulate_viral_rebound_time(new; breakpoint)
			rebound_times[ii,jj] = t_rebound
		end
	end
	return rebound_times
end







## 
# --------------------------------
# Statistical analysis
# --------------------------------



###

##
function concordance_vector(
	t_rebounds, # unadjusted simulations n_patients list of rebound times
	simulation_matrix # (sims vs. trebounds)
	)
	(_, sims) = size(simulation_matrix)
	out = zeros(sims)
	Threads.@threads for ii in 1:sims
		out[ii] = length(t_rebounds) * Tomoko.observed_discrepancy(t_rebounds,simulation_matrix[:,ii])
	end
	#concordance_vector = map(x->length(t_rebounds) * Tomoko.observed_discrepancy(t_rebounds,x),
	#	 eachcol(simulation_matrix))
	return out
end

function generate_statistic_closure(
	mult, # vector of multipliers
	rebsmats... # a list as long as the number of trials
	)
	let null = argmin(ii -> log(mult[ii])^2,  eachindex(mult));
	function stat_fun(t_rebs...)
		total_concord = sum(
			concordance_vector(t_reb,reb_mat) for (t_reb, reb_mat) in Iterators.zip(t_rebs,rebsmats)
			)
		min = argmin(ii -> total_concord[ii], eachindex(total_concord)) 
		(total_concord[null] - total_concord[min],mult[min])
	end
	return stat_fun
	end
end

function generate_concord_vector_closure(
	mult, # vector of multipliers
	rebsmats... # a tuple as long as the number of trials being analyzed together
	# contains distributional information about the 
	)
	let null = argmin(ii -> log(mult[ii])^2,  eachindex(mult));
	function stat_fun(t_rebs...)
		total_concord = sum(
			concordance_vector(t_reb,reb_mat) for (t_reb, reb_mat) in Iterators.zip(t_rebs,rebsmats)
			)
	end
	return stat_fun
	end
end



function discrepancy_stochastic_rebound(trial, ab_list ; 
    data_size=16, min= 1, max=56, # parameters of the data
    samples_per_patient = 10^1,
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0) # * length of θ_vec = generated samples
	t_rebounds = get_observed_rebound_times(trial)
	t_simulated = simulate_trial_rebounds(trial, ab_list; 
		samples_per_patient, selection_multiplier, diversity_multiplier)
	Tomoko.observed_discrepancy(t_rebounds, t_simulated)
end

#----------------------------------------------------------------#
# Rebound probability and
#----------------------------------------------------------------#
using Combinatorics

function rebound_dict(ab_list, n_antibodies; 
		diversity_multiplier = 2.07, # best fit multiplier
		samples_per_patient = 8, 
		quartile_estimator_samples = 20) 
    out = Dict{Array{String,1},Array{Float64,1}}()
    for combo in combinations( ab_list, n_antibodies)
        list = Float64[]
        for ii in 1:quartile_estimator_samples
			p = simulate_trial_rebound_prob_bayes(combo; samples_per_patient, diversity_multiplier)
			push!(list, p)
        end
        out[combo] = list
    end
    return out
end





#----------------------------------------------------------------#
# Make some plots about mutations
#----------------------------------------------------------------#
##
fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_traces.h5", "w")
fid["10-1074"] = trial_traces("10-1074",["10-1074"], samples_per_patient = 100)
fid["3BNC117"] = trial_traces("3BNC",["3BNC117"], samples_per_patient = 100)
fid["combo"] = trial_traces("combo",["10-1074","3BNC117"], samples_per_patient = 100)
close(fid)
##





using CairoMakie

lines(mult,val1 .- mean(val1))
CairoMakie.lines!(mult,val2 .- mean(val2), color = :red)
CairoMakie.lines!(mult,val3 .- mean(val3), color = :blue)
current_figure()

lines(mult,16 .* val1 .+ 14 .* val2 .+ 7 .* val3)

## Collect data at base parameter values 
#----------------------------------------------------------------#
# Statistical analysis required for reservoir stuff
#----------------------------------------------------------------#
mult = 2 .^(-2.5:.05:2.5)


@time rebs1074 = mapreduce(
	m->simulate_trial_rebounds("10-1074", ["10-1074"]; samples_per_patient = 2400, diversity_multiplier=m),
	hcat,  
	mult
	)
@time rebs3BNC = mapreduce(
		m->simulate_trial_rebounds("3BNC", ["3BNC117"]; samples_per_patient = 2400, diversity_multiplier=m),
		hcat,  
		mult
	)
@time rebscombo = mapreduce(
		m->simulate_trial_rebounds("combo", ["3BNC117","10-1074"]; samples_per_patient = 2400, diversity_multiplier=m),
		hcat,  
		mult
	)
fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_diversity.h5", "w")
fid["10-1074"] = rebs1074
fid["3BNC117"] = rebs3BNC
fid["combo"] = rebscombo
close(fid)
##
tempfun = generate_statistic_closure(
	mult, # vector of multipliers
	rebs1074, rebs3BNC, rebscombo
	)

tempfun2 = generate_concord_vector_closure(
		mult, # vector of multipliers
		rebs1074, rebs3BNC, rebscombo
		)
##
mat1 = simulate_trial_rebounds_matrix("10-1074",["10-1074"],samples_per_patient = 2000)
mat21 = simulate_trial_rebounds_matrix("3BNC",["3BNC117"],samples_per_patient = 2000)
mat22 = simulate_trial_rebounds_matrix("3BNC",["3BNC117"],samples_per_patient = 2000)
mat2 = mapreduce(
	(x,y)->vcat(x,y)[sample(1:18,14,replace = false)], hcat,
	eachcol(mat21), eachcol(mat22)
	)
mat3 = simulate_trial_rebounds_matrix("combo",["10-1074","3BNC117"],samples_per_patient = 2000)

@time dist = map(tempfun,
	eachcol(mat1),eachcol(mat2),eachcol(mat3))

(observed, best_fit_multplier) = tempfun(get_observed_rebound_times("10-1074"),get_observed_rebound_times("3BNC"),get_observed_rebound_times("combo"))
concord = tempfun2(get_observed_rebound_times("10-1074"),get_observed_rebound_times("3BNC"),get_observed_rebound_times("combo"))



##
tmat1 = simulate_trial_rebounds_matrix("10-1074",["10-1074"],samples_per_patient = 2000, diversity_multiplier = best_fit_multplier)
tmat21 = simulate_trial_rebounds_matrix("3BNC",["3BNC117"],samples_per_patient = 2000, diversity_multiplier = best_fit_multplier)
tmat22 = simulate_trial_rebounds_matrix("3BNC",["3BNC117"],samples_per_patient = 2000, diversity_multiplier = best_fit_multplier)
tmat2 = mapreduce(
	(x,y)->vcat(x,y)[sample(1:18,14,replace = false)], hcat,
	eachcol(tmat21), eachcol(tmat22)
	)
tmat3 = simulate_trial_rebounds_matrix("combo",["10-1074","3BNC117"],samples_per_patient = 2000, diversity_multiplier = best_fit_multplier)
##

fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_statistics.h5", "w")
fid["multiplier"] = mult
fid["statistic"] = dist
fid["observed_stat"] = observed
fid["best_fit_multplier"] = best_fit_multplier
fid["untransformed/10-1074"] = mat1
fid["untransformed/3BNC117"] = mat2
fid["untransformed/combo"] = mat3
fid["transformed/10-1074"] = tmat1
fid["transformed/3BNC117"] = tmat2
fid["transformed/combo"] = tmat3
fid["concordvec"] = concord
close(fid)


pvalue = 1-ecdf(dist)(observed)

## Collect data at base parameter values 
@time rebs1074 = simulate_trial_rebounds("10-1074", ["10-1074"]; samples_per_patient = 300)
@time rebs3BNC = simulate_trial_rebounds("3BNC", ["3BNC117"]; samples_per_patient = 300)
@time rebscombo = simulate_trial_rebounds("combo", ["3BNC117","10-1074"]; samples_per_patient = 300)
fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_279.h5", "w")
fid["10-1074"] = rebs1074
fid["3BNC117"] = rebs3BNC
fid["combo"] = rebscombo
close(fid)

## mutation effect on rebound probability for the three trials
 rebound_prob_1074 = [simulate_trial_rebound_prob("10-1074",["10-1074"]; samples_per_patient = 10000, no_mutations = false),
 simulate_trial_rebound_prob("10-1074",["10-1074"]; samples_per_patient = 10000, no_mutations = true)];

 rebound_prob_3BNC = [simulate_trial_rebound_prob("3BNC",["3BNC117"]; samples_per_patient = 10000, no_mutations = false),
 simulate_trial_rebound_prob("3BNC",["3BNC117"]; samples_per_patient = 10000, no_mutations = true)];

 rebound_prob_combo = [simulate_trial_rebound_prob("combo",["10-1074","3BNC117"]; samples_per_patient = 20000, no_mutations = false),
 simulate_trial_rebound_prob("combo",["10-1074","3BNC117"]; samples_per_patient = 20000, no_mutations = true)];

 fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_mutvsnomut.h5", "w")
 fid["10-1074"] =  rebound_prob_1074
 fid["3BNC117"] = rebound_prob_3BNC
 fid["combo"] =  rebound_prob_combo
 close(fid)

## Generate the distribution of rebound probabilities

out = rebound_dict(ablist, 2)



 ## Tests and scratchwork








##
vp = initialize_viral_population(.02,ab_list)
@time (t_rebound,trace) =  simulate_viral_rebound_time(vp::ViralPopulation)
##
@btime vp = initialize_viral_population(.1,ab_list)

@btime vp2 = restart_viral_population(vp, [.5,.7,.8])

trace = virus_time_trace(vp, 0:56)

function packbits2(vec::BitVector)
	x::Int64 = 0
	for ii in Iterators.reverse(vec)
		x += ii
		x <<= 1
	end
	x >>=1
	return x
end

function random_int2(frequency_vector::Vector{Float64})
	packbits(rand(length(frequency_vector)) .> frequency_vector)
end

##
lines(0:56, trace[1,:] .+ tracey[2,:], color = :blue)
lines!(0:56, trace[1,:], color = :red)
current_figure()

##