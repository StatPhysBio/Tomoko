using Tomoko
using HDF5
using Statistics


##

#=
In this section we build the results of the analysis of the sites.
This includes the results of the MLE analysis and the Bayesian results.
=#

filepath = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/snpanalysis.h5"
fid = h5open(filepath, "r+")
dictovec(a::Dict) = [a["$i"] for i in 1:length(a)] #auxiliary function for turning dictionaries into Vector
for ab in fid
    global bayesout = []
    for site in ab
        mut = reverse(read(site["mut"]))
        ls = [ [
            Tomoko.LikelihoodSample(mut..., read(day["theta"]), dictovec(read(day["kwt_mut"]))... ; n_samp = 10000) 
            for day in patient if (sum(dictovec(read(day["kwt_mut"]))) > 0) & (read(day["theta"]) > 0)]
        for patient in site["snpdata"] ]
            for avg in [0, 1]
                theta_list = Tomoko.baysian_pop_rs(ls; samples = 2*10^3, burn_in = 10^2, avg = avg)
                site["bayes_avg$(round(avg,digits=1))"] = theta_list
                println("$avg")
            end
    end
end
close(fid)





filepath = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/snpanalysis.h5"
fid = h5open(filepath, "r+")
dictovec(a::Dict) = [a["$i"] for i in 1:length(a)] #auxiliary function for turning dictionaries into Vector
for ab in fid
    for site in ab
        mut = reverse(read(site["mut"]))
        ls = [ [
            Tomoko.LikelihoodSample(mut..., read(day["theta"]), dictovec(read(day["kwt_mut"]))... ; n_samp = 10000) 
            for day in patient if (sum(dictovec(read(day["kwt_mut"]))) > 0) & (read(day["theta"]) > 0)]
        for patient in site["snpdata"] ]
        (rs,g,h) = Tomoko.rs_newton_search(ls)
        site["rsel"] = [rs ,abs(h)^(-1/2)]
        (rs,g,h) = Tomoko.rs_newton_search(ls; average=false)
        site["rsel_ind"] = [rs ,abs(h)^(-1/2)]
    end
end
close(fid)
##
filepath = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/snpanalysis.h5"
fid = h5open(filepath, "r")
textout = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/texfile.tex"
open(textout, "w") do io
    write(io, " \\begin{tabular}{llllcccccc} \n ")
    for ab in keys(fid)
        write(io, "\\hline \\ \\ $ab esc site & wt & mut & \$ μ_{wt} \$ & \$ μ_{mut} \$ & \$ Δ/μ_{ts} \$  & \$ σ_{Δ/μ}\$ & \$ Δ'/μ_{ts} \$  & \$ σ_{Δ'/μ}  \\\\ \n \\hline \n ")
        for site in fid[ab]
            loc = read(site["loc"])
            start = string(read(site["wt_mut_aa"])["1"]...)
            mutat  = string(read(site["wt_mut_aa"])["2"]...)
            mutwt = reverse(read(site["mut"]))[1]
            mutmut = reverse(read(site["mut"]))[2]
            (r,rσ) = read(site["rsel"])
            (r1,rσ1) = read(site["rsel_ind"])
            write(io, " $loc & $start & $mutat & $(round(mutwt,digits = 2)) & $(round(mutmut,digits = 2)) &  $(round(r,digits=1)) & $(round(rσ, digits = 1)) &  $(round(r1,digits=1)) & $(round(rσ1, digits = 1))  \\\\ \n")
        end
    end
    write(io, "\\end{tabular} ")
end;
close(fid)
#$round(mutwt,digits = 2) & $round(mut[2],digits = 2) &  $round(r,digits=1)  $round(rσ, digits = 1) 

