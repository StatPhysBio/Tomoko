using Colors, HDF5, Distributions, DataStructures

path = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/Julia results/"
##
run(`rm $(path)colorlist.h5`)
##

function base_contrast_colors(n,off)
    sort([HSL(mod(MathConstants.φ*ii/2*360+off,360), .7+.3*cos(ii/(n+1)*pi), asin((ii/n)^3)/π +.4) for ii in (1-n):2:n],by = x->x.h)
end


fid = h5open("$(path)colorlist.h5", "cw")

function color_to_tuple(c)
    x = RGB(c)
    convert.(Float64,Float64[x.r, x.g, x.b])
end

aacodons = Dict([["A", ["GCT", "GCC", "GCA", "GCG"]], ["R", ["CGT", "CGC", "CGA", 
   "CGG", "AGA", "AGG"]], ["N", ["AAT", "AAC"]], ["D", ["GAT", 
   "GAC"]], ["C", ["TGT", "TGC"]], ["Q", ["CAA", 
   "CAG"]], ["E", ["GAA", "GAG"]], ["G", ["GGT", "GGC", "GGA", 
   "GGG"]], ["H", ["CAT", "CAC"]], ["I", ["ATT", "ATC", 
   "ATA"]], ["L", ["TTA", "TTG", "CTT", "CTC", "CTA", 
   "CTG"]], ["K", ["AAA", "AAG"]], ["M", ["ATG"]], ["F", ["TTT", 
   "TTC"]], ["P", ["CCT", "CCC", "CCA", "CCG"]], ["S", ["TCT", "TCC", 
   "TCA", "TCG", "AGT", "AGC"]], ["T", ["ACT", "ACC", "ACA", 
   "ACG"]], ["W", ["TGG"]], ["Y", ["TAT", "TAC"]], ["V", ["GTT", 
   "GTC", "GTA", "GTG"]], ["*", ["TAA", "TGA", "TAG"]]])

aatype = Dict(
[["aliphatic",[string(ii) for ii in "AILMV"]],
["aromatic",[string(ii) for ii in "FWY"]],
["polar",[string(ii) for ii in "NCQST"]],
["acidic",[string(ii) for ii in "DE"]],
["basic",[string(ii) for ii in "RHK"]],
["unique",[string(ii) for ii in "GP"]],
["stop",["*"]]])

aaorder = vcat(map(x->aatype[x],["aliphatic","aromatic","polar","acidic","basic","unique","stop"])...)
4^3 == length(vcat(map(x->aacodons[x],aaorder)...)) # covers all codons.
aacolors = reverse(sort(deleteat!(sort(distinguishable_colors(23, transform = deuteranopic), by = x->getfield(HSL(x),:l)),[5,23]), 
    by = x->getfield(HSL(x),:h))[vcat([1,21],collect(2:20))])
aacolorscheme = Dict( aa => color_to_tuple(col) for (aa,col) in zip(aaorder,aacolors))

fid["order"] = aaorder
for key in keys(aacolorscheme)
    fid["colorscheme/$key"] = aacolorscheme[key]
end

for key in keys(aatype)
    fid["types/$key"] = aatype[key]
end

regions = OrderedDict([
    ["V2",Dict{String,Any}("start"=> 158,"end"=>196)],
    ["Loop D",Dict{String,Any}("start"=>275,"end"=>283)],
    ["V3",Dict{String,Any}("start"=>296,"end"=>337)],
    ["V5 & beta23-24",Dict{String,Any}("start"=>451,"end"=>471)],
    ["MPER",Dict{String,Any}("start"=>670,"end"=>699)],
    ["Interface",Dict{String,Any}("start"=>512,"end"=>639)],
    ["CD4 bs",Dict{String,Any}("start"=>364,"end"=>375)]
    ]
)


for (reg, col) in zip(keys(regions),
        base_contrast_colors(length(regions),-95))
    global regions[reg]["color"] = color_to_tuple(col)
end


for key in keys(regions)
    dict = regions[key]
    fid["regions/$key/color"] = dict["color"]
    fid["regions/$key/start"] = dict["start"]
    fid["regions/$key/end"] = dict["end"]
end

close(fid)

##
