export PopRates, PopState, FitVector
export initialize_pop, random_wright_sample, mean_fit

## difinition of the fitness structure

struct FitVector # Custom sparse array type for fitness vector stuff
    β0::Float64
    sind::Array{Int64,1} # sites under selection
    nind::Array{Int64,1} # neutral sites
    svals::Array{Float64,1} # Δ fitness if mutant
    βwt::Float64 # β0 + sum(β1), the wildtype fitness (zero bitvector)
    βmax::Float64 # β0 + sum(abs.(β1)) the maximum fitness
end

function FitVector(β0::Float64, β1::Array)
    # vector of fitness differences between wildtype and mutant
    FitVector(
        β0,
        findall(!iszero, β1),
        findall(iszero, β1),
        -β1[findall(!iszero, β1)],
        β0 + sum(β1)/2,
        β0 + sum(abs.(β1))/2
    )
end


function fitness(f::FitVector, x::BitArray)
    s = 0
    for (ii,val) in zip(f.sind,f.svals)
        if x[ii]
            s += val
        end
    end
    return s + f.βwt
end


# the fitness vector is stored in cachey sparsy form.
# unfortunately we may want to access it in a more conventional way
function β1_vector(f::FitVector)
    # vector of fitness differences f(wildtype) - f(mutant)
    a = zeros(Float64, maximum([maximum([0; f.sind]),maximum([0; f.nind])]))
    for (ii,val) in zip(f.sind,f.svals)
        a[ii] = -val
    end
    return(a)
end


mutable struct PopState{T}        # parameters defining a population state
    individuals::Vector{BitArray{1}}
    fitness::Vector{T}         # tendency toward excess birth
    offset::T # which fitness is currently the breakeven point
    size::Int
    time::T
end


# Parameters defining the simulaiton dynamics
struct PopRates        # parameters defining a population
    loci::Int64        #    = 2^8        number of loci on the genome
    κ::Float64         #    = 1000.        avg pop size
    λ::Float64        #    = 4 (population birth rate per individual (inverse generation time))
    ν::Float64         #    = 1.0e-3    mutation rate per site (per generation)
    χ::Float64        #    = 0.0        outcrossing rate (between 0 and one), rate of recombination events
    ρ::Float64         #    = 1.0e-2    how tightly two crossed genomes get wound together (# of crosses per nucleotide)
    β0::Float64     #    = 1.0        base (competetive) fitness      
    βf::FitVector    #    = FitVector(β0,zeros(loci)) encodes fitness function efficiently
    f::Function     #    = (x->fitness(βf,x))    # mapping from genotype to fitness (increased birth + competion)
    # if you want epistasis, you can code it in by hand using this function.
    function PopRates(loci,κ,λ,ν,χ,ρ,β0,β1)
        # we make βf and f only accessible through inner constructors
        # this enforces their correctness
        if length(β1) == loci
            βf = FitVector(β0,β1)
            f = (x->fitness(βf,x))
            return new(loci,κ,λ,ν,χ,ρ,β0,βf,f)
        else
            error("Selection coefficients wrong length")
        end
    end
end

# default constructor method
function PopRates(;
    loci = 2^8, κ = 1000, λ = 4, ν = 1.0e-3, χ = 0.0, ρ = 1.0e-2, β0 = 1.0,
    β1 = zeros(loci))
    PopRates(loci,κ,λ,ν,χ,ρ,β0,β1)
end

# updating method, pass an old PopRates, and some keywords for parameter changes
function PopRates(par::PopRates; 
    loci = par.loci, κ = par.κ, λ = par.λ, ν = par.ν, χ = par.χ, ρ = par.ρ, β0 = par.β0,
    β1 = β1_vector(par.βf))
    PopRates(loci,κ,λ,ν,χ,ρ,β0,β1)
end


function mean_fit(par::PopRates)
    par.β0 + mean( β1_vector(par.βf)' *(1 .- 2 .* wright_sample(par)) for ii in 1:100)/2
end

## Equilibrium Sampling
# We assume that the population is very near the maximum


function pop_ne(par::PopRates)
    par.κ/par.λ
end

function pop_θ(par::PopRates)
    2*pop_ne(par)*(par.ν*par.λ)/2
end

"""
Sample from Wright equilibrium
"""
function random_wright_sample(θ_wt, θ_mut, σ)
    k = rejection_sample_k(θ_wt, θ_mut, σ)
    rand(Beta(θ_wt + k, θ_mut)) # frequency
end

function wright_sample(par::PopRates)
    freq = zeros(par.loci)
    θ_ref = pop_θ(par)
    ne = pop_ne(par)
    # set the selected frequencies
    for (ii, val) in zip(par.βf.sind, par.βf.svals)
        if val < 0 # if mutant has a disadvantage
            freq[ii] = 1 -random_wright_sample(θ_ref, θ_ref, 2*ne*abs(val))
        else # if the mutant has an advantage
            freq[ii] = random_wright_sample(θ_ref, θ_ref, 2*ne*abs(val))
        end
    end
    for ii in par.βf.nind # for the neutral indices, sample with σ = 0
        freq[ii] = random_wright_sample(θ_ref, θ_ref, 0)
    end
    return freq
end

"""
Generate individuals with site entries equal to 1 according to binomial
frequencies 
"""
function initialize_pop(frequencies::Vector, pop_size::Int, par::PopRates; time = 0.0)
    genomes = mapreduce(p->rand(Bernoulli(p), pop_size),hcat,frequencies)
    individuals = [BitArray(genomes[ii,:]) for ii in 1:pop_size]
    fitness = par.f.(individuals)
    offset = sum(fitness)/par.κ
    pop = PopState(
        individuals,
        fitness,
        offset,
        pop_size,
        time
        )
end

# asumption of neutrality
function initialize_pop(par::PopRates; pop_size = ceil(Int, par.κ))
    initialize_pop(wright_sample(par), pop_size, par) # draw each site from Wright eq distribution
end

# simpler specification for constant frequencies
# useful for indicating pure wildtype
function initialize_pop(freq::Real, par::PopRates; pop_size = ceil(Int, par.κ))
initialize_pop([freq for ii in 1:par.loci], pop_size, par)
end

"""
Recalculate the caches
"""
function renew_fitness(pop::PopState, par::PopRates)
    pop.fitness = par.f.(pop.individuals)
    pop.offset = sum(pop.fitness)/par.κ
    pop.size = length(pop.individuals)
end

