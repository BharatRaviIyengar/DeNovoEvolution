using AbstractTrees
using StatsBase

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

function phyloprob(plist::Dict,anclen::Int,transmat::Matrix{Float64})
    nlist = popnodes[[x.name in keys(plist) for x in popnodes]]
    mrca = LCA(nlist);
    anclist = [vcat(x,lastancestors(x,mrca)...) for x in nlist]
    tpvals = 1
    countedans = []
    for i in eachindex(nlist)
        tblen = sum([x.age for x in anclist[i]])*26
        flen = plist[nlist[i].name]
        if !(isempty(countedans))
            repnodes = (anclist[i])[anclist[i] .∈ Ref(countedans)]
            tblen -= sum([x.age for x in repnodes])*26
        end
        countedans=union(countedans,anclist[i]);
        tpvals*=(transmat^tblen)[anclen-3,flen-3];
    end
    return tpvals
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
            tpvals = [phyloprob(plist,anclen,tmatX) for anclen in luniq]
            mplen = luniq[tpvals .== maximum(tpvals)]
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

datax =  readdlm(outfile)[:,[6,8]];
tmdict = Dict(x => tmatX^x for x in unique(datax[:,2]));

estprobs = [tmdict[x[2]][x[1],x[1]]/sum(tmdict[x[2]][x[1],:]) for x in eachrow(datax)];

uprobrng = rand(size(datax)[1],100000);

estnum =[sum(x.>=estprobs) for x in eachcol(uprobrng)]
