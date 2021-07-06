using Tomoko
using HDF5
using Statistics
using Combinatorics

trials = map(x->"/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialanalysis$x.h5", ["10-1074","3BNC","combo"])
snpanalysis = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/snpanalysis.h5"
## Generate theta vector
begin
global thetas = []
for path in trials
h5open(path, "r") do fid
    for patient in fid["genetic"]
        for day in patient
            if (read(day["day"]) < 2) & haskey(day, "theta")
                push!(thetas,read(day["theta"]))
            end
        end
    end
end
end
end

## retrieve antibody profile list (ablist)
# and dictionary of uncertainties (abrσ)

begin
global ablist = Dict()
global abrσ = Dict()
h5open(snpanalysis, "r") do fid
        abs = keys(fid)
        for ab in abs
            sites = Array{Float64,1}[]
            siteσ = Float64[]
            for site in fid[ab]
                rs = read(site["rsel"])[1]
                mut = reverse(read(site["mut"]))
                rσ = exp((log(read(site["rsel"])[2]) + log(read(site["rsel_ind"])[2]))/2)
                if rs > 0
                    push!(sites,vcat(mut...,rs))
                    push!(siteσ,rσ)
                end
            end
            if ab != "101074"
            ablist[ab] = sites
            abrσ[ab] = siteσ
            end
        end
    end
end

##

function rebound_dict(thetas, ablist, n_antibodies, samples) 
    out = Dict{Array{String,1},Float64}()
    for combo in combinations( collect(keys(ablist)), n_antibodies)
        val = 0
        for θ in thetas
            sample_array = [map(x-> Tomoko.WrightSampler((θ .* x)...), ablist[ab]) for ab in combo]
            val += mean( 56 < -3 * log(
                prod( 1 - prod(Tomoko.sample!.(arr)) for arr in sample_array)
                ) for _ in 1:samples )
        end
        out[combo] = val/length(thetas)
    end
    yy = sort(collect(zip(values(out),keys(out))), by = x->-x[1])
    vals = map(x->x[1],yy)
    combos = map(x->x[2],yy)
    return(vals,combos)
end

function random_sampler_array(θ, x, r_σ)
    Tomoko.WrightSampler((θ .* (x .* exp.([0, 0, (r_σ/x[3]) *randn()])))...)
end

function rebound_dict_σ(thetas, ablist, abrσ, n_antibodies, samples) 
    out = Dict{Array{String,1},Array{Float64,1}}()
    for combo in combinations( collect(keys(ablist)), n_antibodies)
        list = Float64[]
        for ii in 1:100
        val = 0
        for θ in thetas
            sample_array = [map( (x,r_σ)-> random_sampler_array(θ, x, r_σ), ablist[ab], abrσ[ab]) for ab in combo]
            val += mean( 56 < -3 * log(
                prod( 1 - prod(Tomoko.sample!.(arr)) for arr in sample_array)
                ) for _ in 1:samples )
        end
        push!(list, val/length(thetas))
        end
        out[combo] = list
    end
    return out
end

path = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/Julia results/"
cd(path)
for n in 1:3
    uncertainty = rebound_dict_σ(thetas, ablist, abrσ, n, 10^3) 
    (vals, combos) = rebound_dict(thetas, ablist, n, 10^3) 
    h5write("rebounds.h5", "mean/ab$n/vals",vals)
    h5write("rebounds.h5", "mean/ab$n/dist" , hcat(map(x->uncertainty[x],combos)...))
    h5write("rebounds.h5", "mean/ab$n/combos",map(x->string((x.*" ")...)[1:end-1],combos))
end

