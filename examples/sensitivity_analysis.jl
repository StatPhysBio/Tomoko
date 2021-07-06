
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




trialpath(x) = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialanalysis$x.h5"
snpanalysis = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/snpanalysis.h5"

#for trial in ["10-1074","3BNC","combo"])

begin
    global thetas = []
    global escapefreq = []
    let path = trialpath("10-1074")
    h5open(path, "r") do fid
        for patient in fid["genetic"]
            for day in patient
                if (read(day["day"]) < 2) & haskey(day, "theta")
                    push!(thetas,read(day["theta"]))
                end
            end
        end
        for patient in fid["viraemic/patients"]
            push!(escapefreq,read(patient["x"]))
        end
    end
    end # end let
    abprofile = h5open(snpanalysis, "r") do fid
        ab = fid["10-1074"]
        return [ vcat(reverse(read(site["mut"]))...,read(site["rsel"])[1]) for site in ab]
    end
end