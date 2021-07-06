# Tomoko

Pop-gen simulator for Julia.  Tomoko is a forward-time, individual-based, bi-allelic haploid simulator suitable for testing high-dimensional models of fitness, linkage, and recombination. It is named in honor of Tomoko Ohta, pioneer of the nearly-neutral theory in molecular genetics.  Major inspiration for the UI is from the `high_d` simulator in [FFPopSim](http://webdav.tuebingen.mpg.de/ffpopsim/). 

Also included are several tools associated with Wright equilibrium including sampling and estimation of population and site parameters associated with the equilibrium distribution.

Tomoko Ohta has said,

>Genetics is now at a very interesting stage. There are so many interesting questions unanswered and so many ways to test to find answers. Intuition is very important in addressing questions. Nurture your own sensibility and pursue your research and work with confidence.

# Parameters
Genotypes (stored as `BitVector`'s) and individuals are in a 1-1 mapping. We do not keep track of the number of clones as is done in more efficient low-mutation, low-recombination simulation schemes. The `PopState` struct has a field `.individuals` which is just the vector of the `BitVector` genomes, of length equal to the number of individuals.

The dynamics of a population are determined by a type `PopRates`.  `PopRates` has the following default definitions

```julia
struct PopRates		# parameters defining a population
	loci = 2^8      # Number of loci
	κ = 1000.       # avg pop size
	β0 = 1.0        # base (competetive) fitness  
	β1 = zeros(loci)# fitness difference between mutant and wildtype at each site.
	λ = 4         	# total birth death noise
	ν = 1.0e-3      # mutation rate per site per birth μ = ν*λ/2 at capacit, on average with minor fitness effects.
	χ = 0.0         # outcrossing rate (between 0 and one), rate of recombination events
	ρ  = 1.0e-2     # how tightly two crossed genomes get wound together (crosses per nucleotide)
	f::Function     # fitness function of genotype. Defaults to sparse version of β0 + sum(β1 .*  x) 
end
```

`β1` is an array of fitness differences between wildtype (locus = 0) and the mutant (locus = 1).  These fitnesses are implemented in a symmertric (-1,1) way by default with no epistasis, the genome bit at position `i` `σ_i = [0,1]` determines fitness `f = Σ_i (2σ_i-1) β_i`.  Epistasis can be defined by setting `f` to be an arbitrary fitness function when constructing a `PopRates` object

All individuals are assumed to have the same activity `λ` which is the sum of the birth rates and death rates. The birth rates `b = (λ + f - φ)/2`, and the death rates are `b = (λ - f + φ)/2`, so `b + d = λ`. `φ` is something like the chemical potential, equal to `φ = (Σ_i f_i)/κ ` it sets the competetive pressure and the mean growth rate of the population.

In our Gillespie scheme, we can sample an individual at random and then use a Bernoulli trial to determine whether it has experienced a birth or death event. This makes our system very close to a Moran model. The equivalent Wright-Fisher (WF) effective population size is `N_e = κ/(2*λ)`, so our WF-generation time is `1/(2λ)`.

The efficiciency in this scheme relative to WF for computing competition is something like the computational efficiency of using the grand-cannonical (defined with temperature, chemical potential) vs. cannonical ensembles (defined with temperature and particle number) in statistical mechanics. Our individuals are conditionally independent, and competition is only mediated through the fitness offset `φ`. This makes the math (specifically the parent sampling in our case) much easier  because we do not sample from the multinomial vector with weights `w_i = exp(f_i)/(Σ_i exp(f_i)`, a birth event is O(1) instead of an O(n)) operation (weighted sampling).

# Bulk parameters and simulations

In the bialleleic distribution we have the stationary state

```julia
p(x) = exp(σ * x) * x^θ_wt * (1-x)^θ_mut / (x * (1-x)) / z_eq(θ_wt, θ_mut, σ)
```
where `z_eq(θ_wt, θ_mut, σ) = ∫dx exp(σ * x) * x^θ_wt * (1-x)^θ_mut / (x * (1-x))` is the normalizing constant (and moment generating function) for the stationary distribution.

Diversity parameters `θ_mut = 2 κ μ_mut / λ` and `θ_wt = 2 κ μ_wt / λ`, where these are the forward and backward mutation flux per individual. `μ_mut` is the result of imperfect copies at birth, with the fidelity being given by `ν`. With the current form of simulations, detailed balance holds and `θ_wt = θ_mut = θ_ref`.  Therefore in terms of simulation parameters, ignoring variations in fitness (at capacity `f-φ ~ 0`), the birth rate is `λ/2`, the mutation rate in the bulk sense is given by `μ_mut = ν * λ/2`  and `θ = 2 Ne * μ = n * ν`.

The selection exponent is `σ = 2 * N_e * (f1-f0). `Typically what we can most easily measure from genetic statistics is the selection / mutation ratio`r_s = (f1-f0)/μ = (f1-f0)/(ν * λ/2)` so that `σ = θ * r_s`.

Although the `PopRates` object defines the parameters required to actually run simulations (`μ,λ,β0,β1`) they have to be reinterpreted in terms of biological observables . There is an equivalence class of simulations with the same biologically observable parameters defined by 1.) changes in temporal scale (λ). and 2.) changes in the discretization scale defined by the number of individuals `κ` since at capacity `κ = n`.


# Running simulations

To run a simulation, you 
1. define an instance of a `PopRates` object, by specifying where the fields are different from the default values.
2. `run_sim(par,timesteps)` timesteps tells you both how long to run the simulation and at which time-steps to gather statistics.

```julia
par = PopRates(χ=.2, ρ=0.1, μ=10^-4, loci = 512)
df = run_sim(par, 1:5:10000)
```

`run_sim` initializes a population by drawing frequencies from the Wright equlibrium at each locus.  Then the population is propagated forward in time using an exact Gillespie simulatior while statistics at the specified timepoints. These statistics and time point are stored as a dataframe.

The simulation is in the form of an individual-based chemical reaction model. The time between birth/death events is exponentially distributed. This leads to a stochastic Lotka-Voltera equation for the mean frequency of a particular trait in the absence of linkage.

# The default statistics
We collect the following default statistics as columns in a DataFrame
`[:time, :pop_size, :freq, :D, :mean_fit, :var_fit]`

* `time` The time at which the statistics were collected. Interconvertable with generations in the Wright-Fisher sense
* `pop_size` Number of individuals at a given time
* `freq` vector of frequency of mutants (i.e. ones) at a given loci
* `θ` The genotypic diversity of the population. (technically defined as the maximum likelihood estimator for the dispresion of the beta-binomial distribution of neutral site allellic counts.)
* `mean_fit` The mean fitness of the population
* `var_fit` The variance of the fitness of the population

Plotting, especially with `Gadfly.jl` and `StatPlots.jl` is well integrated with data frames.  For example, using the output of the simulation above, you can plot the `pop_size` over time

```julia
using Gadfly
plot(df, x = :time, y = :pop_size, Geom.line)
```

# Defining more complicated simulations.

In the simulation defined, we start with a populaiton in pseudo-equilibrium and observe the (mild) effects of linkage and drift when the simulation runs. A pseudo-equilibrium starting point is useful if we are interested in the population stationary state because we have to wait less time for things to settle down. However, you are probably running exact simulations because you are interested in what happens *far* from stationary state. This section is about how to spice things up a bit.

## What is inside `run_sim`

Let's just look at what `run_sim` is made of. It's nothing much more than a for loop. Here's the main definition

```julia
function run_sim(pop::Popstate, par::PopRates, timepoints;
	var_sites = [], # vector of indices that are flipping
	flip_prob = 0.1) # probability that they flip at every time point.
    df = initialize_stats() # empty data frame for holding the stats
    for t in timepoints
        run_until!(pop,par,t) # evolve the poulation state until pop.time > t.
        record_stats!(df, pop, par)  # record the statistics in the data_frame
        par = selection_flip(par; var_sites = var_sites, prob = flip_prob) # possibly stochastically change the environment.
    end
    return df # return the df, now heavy with juicy statistics for you to plot or analyze
end

# definition when pop is not provided: linkage-free equilibrium initial condition.
function run_sim(par::PopRates, timepoints; kwdargs...) # population loci set to wright equilibrium, 
	pop = initialize_pop(par)
	run_sim(pop, par, timepoints; kwdargs...)
end

```
Complications can fall into three categories
1. Change start population
1. Vary environment or add external events
1. Calculate more statistics.


## Changing the start population
To make things more intersting, we can change the initial population to be more out-of equilibrium. For instance, we might define the population to be entirely made up of wildtype clones

```julia
pop = initialize_pop(0.0,par) # the first argument specifies the mutant frequency
```

Or we can grow out the population out from just 10  wt-individuals:

```julia
par = PopRates(χ=.2, ρ=0.1, μ=10^-4)
pop = initialize_pop(0.0, par; pop_size = 10) 
df = run_sim(pop, par, 1:5:10000)
```
We can also define a population stochastically with a full vector of frequencies of length equal to the number of loci to be sampled binomially.

## Varying the evironment, bottleneck events.
As seen above, `run_sim` has the ability to include environmental variation, consisting of selection sign flips on some number of active sites through the defined `active_sites` vector.  This gets used by `selection_flip` which returns a new ParRates object.

One can replace where `selection_flip` occurs in the `for` loop with more complicated changes to the evolutionary parameters or with functions on the population, like bottleneck events that remove individuals based on a particular sequence of loci.

Note: `par::PopRates` is an immutable with internal constructors and isn't be mutated on the fly. Instead a new instance is defined and overwrites the local variable `par`. On the other hand,`pop::PopState` is a mutable and can be mutated at will.

## Adding more statistics
The statistics are functions of the population and parameters that live inside the Tomoko module.  This is how `record_stats` identifies the statistic name with the statistic function.  This means that to define more functions from the REPL you have to `eval` them into the Tomoko context and update the global variable `pop_stats` so that `record_stats` knows you want to keep track of a new variable.

```
julia> Tomoko.eval(:(
function new_stat(pop::PopState, par::PopRates)
    do_stuff(pop::PopState, par::PopRates, sites = sites_of_interest)
end
))

julia> push!(Tomoko.pop_stats,:new_stat)
```

 It's cool that Julia lets you do this, but there's a tradeoff, and that is your session becomes dependent on your `eval`-history. From a functional programming persepctive, this is a bad idea. This design was chosen to keep the statistic-gathering machinery as global variables to avoid having to pass yet another argument to the simulation functions and to make sure that the name and function are intrinsincally linked.

All in all, the advantages of the current design, memory and siplicity of commands are outweighed by flexibility of design and fineness of control offered by alternatives. On the memory footprint for reasonable (read: useful) `PopState`'s ended up being the same order of magnitude as that for the sufficient statistics: 1000 or so (individuals), with 512-bit genomes is not much worse than 512 and change statistics, made up of 64 length floats.  In the future, the `PopState` post-processing statistics construction machinerery will be seperated from the simulation machinery.

# Extending and contributing
In the end it's impossible to design a user interface that can do everything from the REPl.  To run the experiment you need to run to answer your scientific questions, you will probably have to look at the source and see what's there, dev and modify the package to suit your needs.

By understanding what methods are available out of the box, you can get a feel for how the machinery works and how to extend it. Help Tomoko.jl evolve with PR's and feature requests!