using Tomoko
using Random
using Test
using Distributions

p = PopRates()
x = bitrand(p.loci)
df = runsim(p,1:1:1000)
sample_betabin(α,β,n,N) = mapreduce(x-> [x,n-x],hcat,rand(BetaBinomial(n,α,β),N))
αtest = 0.05
dat = [sample_betabin(αtest,αtest,100,200) for ii in 1:20]

@testset "Tomoko.jl" begin
    0.0 < p.f(x)
    (mean(df.D) - Tomoko.pop_D(p))/Tomoko.pop_D(p) < 0.05 # check that D is near the equilibrium value
    abs((mean(df.pop_size) - p.κ)/p.κ) < 0.05 # check that the pop_size approaches the eq value
    abs((median(Tomoko.neutral_D.(dat)) - 2 αtest)/(2*αtest)) < 0.1 # check that minka algorithm is correct
end