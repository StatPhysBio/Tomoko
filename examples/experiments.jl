

using Revise
using Distributed
addprocs(8)
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using Tomoko
using Tomoko
using CSV
using Printf
using DataFrames
using Statistics
using HDF5


path = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/Julia results/"
cd(path)
fid = h5open("linkage_sims", "w")

function save_pop_fit(path, result, par; name_fields = [:κ,:ν,:χ], p = 0)
    s = ["$(name)=$(@sprintf("%.0E", getfield(par,name)))_" for name in name_fields]
    CSV.write(string(path,s...,"pflip=$(@sprintf("%.0E", p))",".csv"),result)
end

## Observables 
# θ = 2*μ*Ne  = 2*(ν *λ / 2) κ/λ = ν κ
# rs = Δ/μ  = Δ/(ν *λ / 2)  =>  Δ = (ν *λ / 2) rs

global loci = 2^9
global lambda = 4
global kappa = 1000
for x in [0,0.1], p in [0,0.01], θ_ in 0.02
    β1 = zeros(loci)
    nu = θ_ / kappa # the first number sets the diversity 
    fixed_sites = 2^3:2^3:2^9 # 64 fixed sites running from Δ/μ = 8 to 512 in increments of 8
    for ii in fixed_sites
        β1[ii] += ii * (nu * lambda) /2
    end
    variable_sites = 2^4-1:2^4:2^9 # 32 variable sites
    for ii in variable_sites
        β1[ii] += .02 # 2% fitness difference (max - mean < 4 to avoid blow up)
    end
    true_fit =β1[fixed_sites]/((nu * lambda) /2)
    result = DataFrame(f = true_fit)
    par = PopRates(
        κ= kappa,
        χ= x, 
        ρ= 0.1, 
        ν= nu , # first number determines D
        loci = loci, 
        β1 = β1)
    for ii in 1:80
        dglist = pmap(x->run_sim(par, 0:10:2000 ; var_sites = variable_sites, flip_prob = p/length(variable_sites)), 1:11)
        fitness = pmap(ii->estimate_rs(dglist, locus = ii, timepoints = 1000:100:2000) , fixed_sites)
        result[Symbol("run",ii)]=fitness
    end
    #save_pop_fit(path, result, par, p = p)
end


##

##  Sensitivity analysis 

# test

b = Bins(1 .+ randn(10))

function norm_likelihood(mean, b::Bins; samples = 10^5)
    b.count .= 0 
    for _ in 1:samples
        x = randn() + mean
        count!(x,b)
    end
    hellinger(b)
end

meanvals = -2:.2:2
likvals = [norm_likelihood(ii,b) for ii in meanvals]

plot(meanvals,likvals)

##

pars1074 = [
        [0.7639963491553017, 0.38199817457765084, 505.9257358559],
        [2.5279926983106034, 0.2974109056836004, 302.00565652011676],
        [1.3079957406811853, 0.34121628017770056, 168.57885310136976]]

pars1074 = map(x->x[[2,1,3]], pars1074)

pars3BNC = [
    [1.0659990872888254, 1.0659990872888254, 125.715013338639], 
    [1.1979972618664763, 0.47919890474659055,222.1651518113006], 
    [0.5, 0.5, 126.83837392646784]
]

pars3BNC = map(x->x[[2,1,3]], pars3BNC)

## 10-1074 clinical trial sensitivity_analysis

logmutfreqs = [6.81132, 7.87891, 3.8231, 7.892, 8.15017, 3.69288, 7.74824, 5.44809, 
4.70357, 9.13302, 21.4665, 8.22195, 11.1314, 0, 19.5888, 5.06494]
startpopsize = [0.0192535, 0.0270581, 0.0150585, 0.0176917, 0.00782702, 0.0191018, 
0.0216631, 0.0143935, 0.0138196, 0.00922108, 0.00694056, 0.00710334,
0.0165255, 0.00862304]

bins1074 = Bins(logmutfreqs, min = 0.1, max  = 56/3)
discr = map(x->
    data_discrepancy(bins1074, (10^(x) * 2) .* startpopsize, pars1074),-1:.1:1)
discr2 = map(ii->
    data_discrepancy(bins1074, 2 .* startpopsize, map(x -> x .* [1,1,10^ii],pars1074))
    ,-1:.1:1)
discmap = pmap(
    ((ii,jj),) -> data_discrepancy(bins1074, (10^(ii) * 2) .* startpopsize, map(x -> x .* [1,1,10^jj],pars1074)),
    Iterators.product(-1:.1:1,-1:.1:1))
xvals = map(((ii,jj),) -> ii, Iterators.product(-1:.1:1,-1:.1:1))
yvals = map(((ii,jj),) -> jj, Iterators.product(-1:.1:1,-1:.1:1))
statlist = hellinger_stat_simulations(2 .*startpopsize, pars1074)
observedstat = data_discrepancy(bins1074, 2 .* startpopsize, pars1074)

    path = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/Julia results/"
    cd(path)
    using HDF5
    h5write("test.h5", "grp1074/discmap", discmap)
    h5write("test.h5", "grp1074/xvals", xvals)
    h5write("test.h5", "grp1074/yvals", yvals)
    h5write("test.h5", "grp1074/statlist",  convert(Array{Float64,1}, statlist))
    h5write("test.h5", "grp1074/observedstat", observedstat)


# 3BNC trial sensitivity


logmutfreqs = [0, 2.58509, 0, 3.47146, 0, 3.05345, 0, 0, 5.99143, 20.2276, 5.68782,
20.9857, 7.52069, 8.87487, 6.41905, 9.22242, 20.1416]
startpopsize = [0.0221461, 0.0153592, 0.0124991, 0.0224002, 0.010789, 0.00613182,
0.00745136, 0.0179552, 0.0040067]

bins3BNC = Bins(logmutfreqs, min = 0.1, max  = 56/3)
discr = map(x->
    data_discrepancy(bins3BNC, (10^(x) * 2) .* startpopsize, pars3BNC),-1:.1:1)
discr2 = map(ii->
    data_discrepancy(bins3BNC, 2 .* startpopsize, map(x -> x .* [1,1,10^ii],pars3BNC))
    ,-1:.1:1)
discmap = pmap(
    ((ii,jj),) -> data_discrepancy(bins3BNC, (10^(ii) * 2) .* startpopsize, map(x -> x .* [1,1,10^jj],pars3BNC)),
    Iterators.product(-1:.1:1,-1:.1:1))
xvals = map(((ii,jj),) -> ii, Iterators.product(-1:.1:1,-1:.1:1))
yvals = map(((ii,jj),) -> jj, Iterators.product(-1:.1:1,-1:.1:1))
statlist = hellinger_stat_simulations(2 .*startpopsize, pars3BNC)
observedstat = data_discrepancy(bins3BNC, 2 .* startpopsize, pars3BNC)

    h5write("test.h5", "grp3BNC/discmap", discmap)
    h5write("test.h5", "grp3BNC/xvals", xvals)
    h5write("test.h5", "grp3BNC/yvals", yvals)
    h5write("test.h5", "grp3BNC/statlist",  convert(Array{Float64,1}, statlist))
    h5write("test.h5", "grp3BNC/observedstat", observedstat)


# Combo sensitivity analysis


logmutfreqs = [22.2362, 3.25652, 22.8749, 5.49538, 9.70821, 23.0945, 20.1188]
startpopsize = [0.00316398, 0.0232185, 0.0368048, 0.00878445, 0.00359203, 0.03172,
0.00801984]

binscombo = Bins(logmutfreqs, min = 0.1, max  = 56/3)
discr = map(x->
    data_discrepancy(binscombo, (10^(x) * 2) .* startpopsize, pars1074, pars3BNC),-1:.1:1)
discr2 = map(ii->
    data_discrepancy(binscombo, 2 .* startpopsize, map(x -> x .* [1,1,10^ii],pars1074), map(x -> x .* [1,1,10^ii],pars3BNC))
    ,-1:.1:1)
discmap = pmap(
    ((ii,jj),) -> data_discrepancy(binscombo, (10^(ii) * 2) .* startpopsize, 
    map(x -> x .* [1,1,10^jj], pars1074), 
    map(x -> x .* [1,1,10^jj], pars3BNC)
     ),
    Iterators.product(-1:.1:1,-1:.1:1))
xvals = map(((ii,jj),) -> ii, Iterators.product(-1:.1:1,-1:.1:1))
yvals = map(((ii,jj),) -> jj, Iterators.product(-1:.1:1,-1:.1:1))
statlist = hellinger_stat_simulations(2 .*startpopsize, pars1074, pars3BNC)
observedstat = data_discrepancy(binscombo, 2 .* startpopsize, pars1074, pars3BNC)

    h5write("test.h5", "grpcombo/discmap", discmap)
    h5write("test.h5", "grpcombo/xvals", xvals)
    h5write("test.h5", "grpcombo/yvals", yvals)
    h5write("test.h5", "grpcombo/statlist",  convert(Array{Float64,1}, statlist))
    h5write("test.h5", "grpcombo/observedstat", observedstat)

##
# 95% confidence level = 0.2717
observedstat = data_discrepancy(binscombo, 2 .* startpopsize, pars1074, pars3BNC)
mean(observedstat .< statlist)



using Plots
plot(-1:.1:1 , discr)

plot(-1:.1:1 , discr2)

