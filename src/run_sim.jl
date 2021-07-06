
export run_until!, next_event!

function  gillespie_step(pop::PopState, par::PopRates)
    total = pop.size * par.λ # total rate
    Δt = randexp() / total
    ind = rand(1:pop.size) # choose an individual at random
    if  (2*rand() - 1)*par.λ + pop.fitness[ind] - pop.offset > 0
            process = ind
    else
            process = -ind
    end
    return (process, Δt)
    # sign of process tells us whether it's a birth or death
end

# Mutation

function flipat!(bits::BitArray, inds)
    for ii in inds
        bits[ii] = !bits[ii]
    end
end

"""
Mutates the genome at a rate given by μ
Assumes the action of a markov chain so that 
mutate!(X, 2μ) = mutate!(mutate!(X, μ), μ).
mutate!(Inf) should randomize the genome
"""

function mutate!(genome::BitArray, ν::Real)
    N = length(genome)
    p = -expm1(-2ν)/2
    flips = randsubseq(1:N, p)
    flipat!(genome, flips)
end

"""
Recombines two genomes at a rate given by ν
Assumes the action of a markov chain so that 
mutate!(X, 2ν) = mutate!(mutate!(X, ν), ν).
mutate!(Inf) should randomize the genome

To undertand the code, here we have randomly chosen cuts at [2,5].
From positions 3-5 the codes are swapped as in this diagram.
*--. ,---. ,----->
    /     /
*--' `---' `----->

"""


function recombine!(v1::T, v2::T, ρ::Real) where T <: AbstractArray
    N = length(v1)
    cut = sort(randsubseq(0:(N-1),ρ))
    if isodd(length(cut))
        push!(cut,N)
    end
    for ii in 1:2:length(cut)
         for jj in cut[ii]+1:cut[ii+1]
             (v2[jj], v1[jj]) = (v1[jj], v2[jj])  
         end
    end
end

"""
Add an individual to the population.
Adjust the birth and death rates accordingly
"""
function add_individual!(pop::PopState, child::BitArray{1}, par::PopRates)
    # add to the fitness
    push!(pop.individuals, child)
    push!(pop.fitness, par.f(child))
    pop.offset += last(pop.fitness)/par.κ      # add to competition
    pop.size += 1 
end

"""
Delete an individual to the population.
Adjust the birth and death rates accordingly
"""
function delete_individual!(pop::PopState, index::Int, par::PopRates)
    # add to the death rate
    pop.offset +=  -pop.fitness[index]/par.κ
    deleteat!(pop.fitness, index)
    deleteat!(pop.individuals, index)
    pop.size += -1 
end

function birth_event!(pop::PopState, parent_index::Int, par::PopRates)
    child = copy(pop.individuals[parent_index])
    mutate!(child, par.ν)
    # add to the population
    add_individual!(pop, child, par)
end

death_event!(pop::PopState, index::Int, par::PopRates) = delete_individual!(pop, index, par)


"""
Choose an individual from the population (daddy) to recombine with the individual at
mommy_index.
"""

function recombination_event!(pop::PopState, mommy_index::Int, par::PopRates)
    daddy_index = rand(1:pop.size)
    # update their genomes
    @views recombine!(pop.individuals[mommy_index], pop.individuals[daddy_index], par.ρ)
    # the fitness differences
    Δmfit = par.f(pop.individuals[mommy_index]) - pop.fitness[mommy_index] 
    Δdfit = par.f(pop.individuals[daddy_index]) - pop.fitness[daddy_index]
    pop.fitness[mommy_index] += Δmfit
    pop.fitness[daddy_index] += Δdfit
    # update offset
    pop.offset += (Δmfit+Δdfit)/par.κ
end

function next_event!(pop::PopState, par::PopRates)
    # step forward in time
    (index, Δt) = gillespie_step(pop, par)
    if index > 0
        birth_event!(pop, index, par)
        if rand(Bernoulli(par.χ))
            recombination_event!(pop, index, par)
        end
    else
        death_event!(pop, abs(index), par)
    end
    pop.time += Δt
end

"""
Measure diversity of neutral sites
    data is a matrix of counts ( 2 x n)
"""

function run_until!(pop::PopState, par::PopRates, t) 
    renew_fitness(pop,par)
        # really only neccesary if the parameters have changed
        # slows things down if timesteps are recorded at order generation time
    while pop.time < t
        next_event!(pop, par) 
    end
end




# function run_sim(pop::PopState, par::PopRates, timepoints; 
#     stats = [freq, D, mean_fit, var_fit], 
#     names = [:freq, :D, :fit, :fitvar])
#     rows = []
#     time = 0 
#     for t in timepoints
#         while time < t
#             time += next_event!(pop, par) 
#         end
#         push!(rows,[t;[s(pop,par) for s in stats]])
#     end
#     columns = Vector{Vector}(undef, 0)
#     for (i,r) in enumerate(first(rows))
#         column = Vector{typeof(r)}(undef, length(rows))
#         column .= getindex.(rows, i)
#         push!(columns, column)
#        end
#     df = DataFrame(columns, [:t;names], copycols=false)
# end


##