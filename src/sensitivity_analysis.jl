export Bins, count!, hellinger, data_discrepancy, hellinger_stat_simulations

struct Bins
    Z::Array{Function} # lenght(z) = n List of membership functions
    q::Array{Float64} # lenght(q) = n empirical probabilities associated with each interval
    count::Array{Int64} # simulation counts in each bin
    Bins(Z,q,count) = length(Z) == length(q) ==length(count) ? new(Z,q,count) : error("length mismatch")
end

function Bins(data; min = -Inf, max = +Inf) # censoring cutoffs of continuous measurement and number of measurements outside cutoff
    lows = sum(data .< min)
    highs = sum(data .>= max)
    datapoints = sort(filter(x -> min <= x < max, data)) # get things in order
    midpoints = [(datapoints[ii] + datapoints[ii+1])/2 for ii in 1:length(datapoints)-1] 
    boundaries = cat(-Inf, min, midpoints, max, Inf; dims =1)
    Z = [x -> boundaries[ii] <= x <= boundaries[ii+1] for ii in 1:length(boundaries)-1]
    q = cat(lows, ones(length(datapoints)), highs; dims=1)
    q = q/sum(q)
    count = zeros(Int64,length(q))
    Bins(Z,q,count)
end

function count!(x,b::Bins)
    ind = findfirst(z -> b.Z[z](x) ,eachindex(b.Z))
    b.count[ind] += 1
end

function hellinger(b::Bins)
    total = sum(b.count)
    sum( (sqrt(q) - sqrt(c/total))^2 for (c,q) in zip(b.count,b.q))
end

function observed_discrepancy(reference_data, simulated_data; min = -Inf, max = +Inf)
    bins = Bins(reference_data)
    for x in simulated_data 
        count!(x,bins)
    end
    hellinger(bins)
end


function data_discrepancy(bins::Bins, θ_vec, ab_profiles... ; samples = 10^4) 
    # the data are expected to be minus-log-mutant-frequency objects.
    # The cutoffs are the minimum and maximum data observed.
    b = bins
    b.count .= 0
    for θ in θ_vec
        sampler_array = map(x->WrightSampler(x...),
            hcat([θ .*  profile for profile in ab_profiles]...))
        for _ in 1:samples
            x = -log(escape_freq(sample!.(sampler_array)))
            count!(x, b)
        end
    end
    hellinger(b)
end


function hellinger_stat_simulations(θ_vec, ab_profiles... ; 
    data_size=16, min=0.1, max=56/3, # parameters of the data
    samples = 10^3, # parameters of the simulation
    trials =10^3) # * length of θ_vec = generated samples
    datastream = []
    for θ in θ_vec
        sampler_array = map(x->WrightSampler(x...),
                hcat(θ .* [ ab_profiles...]...)) 
        for ii in 1:trials*data_size
        push!(datastream, -log(escape_freq(sample!.(sampler_array))))
        end
    end 
    out =[]
    while length(datastream) >= data_size 
        data = [];
        for jj in 1:data_size
            ind = rand(1:length(datastream))
            push!(data, splice!(datastream, ind))
        end
        b = Bins(data; min=min,max=max)
        push!(out, data_discrepancy(b, θ_vec, ab_profiles... ; samples = samples))
    end
    return out
end

