using AbstractTrees
using StatsBase
using DelimitedFiles

trimerfreq = readdlm(joinpath(Base.source_dir(),"dmel_intergenic_trimers.txt"), '\t');
nsub, nucsmbt = readdlm("dmel_mutbias.txt",'\t',header = true);
trimerfreq[:,2] = normalize(trimerfreq[:,2]);
stopcodons = ["TAA","TAG","TGA"]
allcodons = kmers(3)
nostopcodons = allcodons[allcodons .∉ Ref(stopcodons)];
noATG = allcodons[allcodons .!="ATG"];


function ATGprobsX(t3)
    ATGprob = nprob3("ATG",t3)
    ATGgain = featuregain3(noATG,["ATG"],t3)
    ATGloss = featuregain3(["ATG"],noATG,t3)/ATGprob
    ATGstay = featurestay3(["ATG"],t3)
    return [ATGprob, ATGgain, ATGloss, ATGstay]
end

function stopprobsX(t3)
    stopprob = sum([nprob3(x,t3) for x in stopcodons])
    stopgain = featuregain3(nostopcodons,stopcodons,t3)
    stoploss = featuregain3(stopcodons,nostopcodons,t3)/stopprob
    stopstay = featurestay3(stopcodons,t3)
    return [stopprob, stopgain, stoploss, stopstay]
end

function tprobs(i,j,ATG,stop)
    nostopstay = 1 - stop[xGain] - stop[xProb]
    noATGstay = 1 - ATG[xGain] - ATG[xProb]
    if i==j    
        return (1-ATG[xLoss])*(1-stop[xLoss])*(1-stop[xGain]/stop[xProb])^i
    elseif i<j
        return (stop[xLoss]*stop[xProb] + ATG[xGain])*nostopstay^(j-i) +(j-i)*stop[xLoss]*stop[xProb]*ATG[xProb]*nostopstay^(j-i-1)
    else
        return stop[xGain]/(1-stop[xProb]) + ATG[xLoss]*noATGstay^(i-j)*ATG[xProb] + ATG[xProb]*(i-j)*stop[xGain]
    end
end

ncodons =[3:900;];
atgvalsX = ATGprobsX(trimerfreq);
stopvalsX = stopprobsX(trimerfreq);
tmatX = [tprobs(i,j,atgvalsX,stopvalsX) for i in ncodons, j in ncodons];

mutable struct datedNode
    name::String
    age::Int64
    children::Vector{datedNode}
    parent::Union{Nothing,datedNode}
end

AbstractTrees.nodevalue(x::datedNode) = x.name*"("*string(x.age)*")"
AbstractTrees.children(x::datedNode) = x.children
AbstractTrees.ParentLinks(::Type{<:datedNode}) = StoredParents()

AbstractTrees.parent(x::datedNode) = x.parent

function newnode(name::String, age::Int64)
    return datedNode(name,age,[],nothing)
end

function addchild(parent::datedNode, child::datedNode)
    child.parent = parent
    push!(parent.children, child)
end

function ancestors(node::datedNode)
    anc = []
    cnode = node
    parent = cnode.parent
    while(!(isnothing(parent)))
        push!(anc,parent)
        cnode = parent 
        parent = cnode.parent
    end
    return anc
end

function lastancestors(node::datedNode,lca::datedNode)
    anc = []
    cnode = node
    parent = cnode.parent
    while(parent != lca)
        push!(anc,parent)
        cnode = parent 
        parent = cnode.parent
    end
    return anc
end

function LCA(x::Vector{datedNode})
    first(intersect(ancestors.(x)...))
end

nodename = (x) -> [y.name for y in x]

# dmelnw = "((GI:9787,(((AK:6400,UM:6400):1180,(DK:6204,SW:6204):1375):638,YE:8217):1570):3056,ZB:12843);"

dmeltree = newnode("Anc",0);
addchild(dmeltree,newnode("Zamb",12843));
SCN = newnode("SCN",1375);
ENE = newnode("ENE",1180);
addchild(ENE,newnode("AK5",6400));
addchild(ENE,newnode("UM",6400));
addchild(SCN,newnode("DK5",6204));
addchild(SCN,newnode("SW5",6204));
NE = newnode("NE",638);
addchild(NE,SCN);
addchild(NE,ENE);

CNE = newnode("CNE",1570);
addchild(CNE,NE);
addchild(CNE, newnode("YE",8217));
EU = newnode("E",3056);
addchild(EU,CNE);
addchild(EU,newnode("GI5",9787));
addchild(dmeltree,EU);

popnodes = collect(Leaves(dmeltree));

function phyloprob(plist::Dict,anclen::Vector{Int},transmat::Matrix{Float64})
    ancidx = anclen .- ncodons[1]
    nlist = popnodes[[x.name in keys(plist) for x in popnodes]];
    mrca = LCA(nlist);
    anclist = [vcat(x,lastancestors(x,mrca)...) for x in nlist]
    countedans = []
    tblen, flidx = [zeros(Int,size(nlist)) for i = 1:2]
    tpvals = ones(size(transmat,1))
    for i in eachindex(nlist)
        tblen[i] = sum([x.age for x in anclist[i]])
        flidx[i] = plist[nlist[i].name] - ncodons[1]
        if !(isempty(countedans))
            repnodes = (anclist[i])[anclist[i] .∈ Ref(countedans)]
            tblen[i] -= sum([x.age for x in repnodes])
        end
        countedans=union(countedans,anclist[i]);
        # tpvals=tpvals.*(transmat^tblen[i]*26)[:,flidx[i]];
    end
    tpvals = .*([(tmatX^(tblen[i]*26))[:,flidx[i]] for i in eachindex(tblen)]...);
    return ancidx[argmax(tpvals[ancidx])]+ncodons[1]
end

file = homedir()*"/Documents/ORF-length-evol/Merged_File_HOM_ORF_wSynteny.txt"

outfile = homedir()*"/Documents/ORF-length-evol/All_summary_wSynteny.txt"

open(outfile,"w") do fout
    for line in eachline(file)
        x = split(line,"\t");
        println(x[1])
        names = split(x[2],",");
        lens = parse.(Int,split(x[3],","));
        luniq = unique(lens)
        local plist = Dict(names[i] => lens[i] for i in eachindex(names))
        if x[6] == "Change"
            mplen = phyloprob(plist,luniq,tmatX)
        else
            mplen = luniq
        end
        medlen = median(lens)
        local ancplist = LCA(popnodes[[x.name in keys(plist) for x in popnodes]])
        divtime = 12843*26
        if !(isempty(ancestors(ancplist)))
            divtime -= (sum([x.age for x in ancestors(ancplist)]) + ancplist.age)*26
        end
        transtat = replace(x[7], "Transcript" => 1, "Not_transcribed" => 0)
        pline = join([x[6],x[2:5]...,mplen[1],medlen,divtime,transtat],"\t")
        println(fout,pline)
        
    end
end

## Maximum Parsimony Inference ##

datax =  readdlm(outfile)[:,[6,8]];
tmdict = Dict(x => tmatX^x for x in unique(datax[:,2]));

estprobs = [tmdict[x[2]][x[1],x[1]]/sum(tmdict[x[2]][x[1],:]) for x in eachrow(datax)];

uprobrng = rand(size(datax)[1],100000);

estnum =[sum(x.>=estprobs) for x in eachcol(uprobrng)]

datay = readdlm(outfile)[:,[1,2,3,6]];
x = findall(datay[:,1] .== "Change");

lengths = [parse.(Int,split(datay[w,3],",")) for w in x];
names = [split(datay[w,2],",") for w in x];
popnames = [k.name for k in popnodes];

lgt2 = [length(x)>2 for x in lengths]

lengths = lengths[lgt2]
names = names[lgt2]
ML_alen = datay[x,4][lgt2];

MP_alen = zeros(Int, size(lengths))
q = 0 
for j in eachindex(lengths)
    cts = countmap(lengths[j]);
    mrep = maximum(values(cts));
    if sum(values(cts) .== mrep) == 1
        MP_alen[j] = maximum(cts).first
    else
        if all(values(cts) .==1)
            MP_alen[j] = 0
            
        else
            q +=1
            mreplens = collect(keys(cts))[values(cts) .== mrep];
            mrepnames = names[j][lengths[j] .∈  Ref(mreplens)]
            aa = LCA(popnodes[[k.name in mrepnames for k in popnodes]]);
            outgrp = aa.children[aa.children .∈  Ref(popnodes)]
                # if isempty(outgrp)
                #     outgrp
                # end
            outlen = names[j] .== outgrp[1].name;
            MP_alen[j] = lengths[j][outlen][1]
        end
    end
end

