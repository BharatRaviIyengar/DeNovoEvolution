using Plots
cd(Base.source_dir())
include("nucleotidefuncts.jl")
stop = ["TAA","TAG","TGA"]

function nostop(codonset)
	codonset[codonset .∉ Ref(stop)]
end

gccontent = 0.493611; # all mRNA #
# gccontent = 0.53942 # all CDS #

p_T = 0.204217;
p_A = 0.256362;
p_G = 0.267809;
p_C = 0.271612;

allcodons = kmers(3)
nostopcodons = nostop(allcodons);

syncodons = Dict(x => gencode[(gencode[:,3].== transnucs(x)) .& (gencode[:,1] .!= x),1] for x in nostopcodons)

PBMEC = 0 .+readdlm("PMBEC.txt")[2:21,2:21]
PBMEC2 = PBMEC .- 0.3;
setindex!.(Ref(PBMEC2), 0.0, 1:20, 1:20);

aanames = Dict(aas[i] => i for i in 1:20)

simcodons = Dict{String,Vector{AbstractString}}()
for x in keys(syncodons)
    y = transnucs(x)
    simaas = aas[PBMEC[:,aanames[y]] .> 0.05]
    simaas = simaas[simaas .!= y]
    simcody = gencode[gencode[:,3] .∈ Ref(simaas),1]
    simcodons[x] = vcat(syncodons[x],simcody)
end

vjoin = (z) -> [join(x) for x in eachrow(z)]
phcp = (x,y) -> permutedims(hcat(collect.(product(x,y))...))

p_synsubs = [featuregain([x],syncodons[x],gccontent) for x in keys(syncodons)]
p_simsubs = [featuregain([x],simcodons[x],gccontent) for x in keys(simcodons)]

# Codon pairs that have a stop codon in frame -2
# Example: 
# TTCACG
#  AGT
#
# Codon pairs that have a stop codon in frame -1
# Example: 
# TCTCAG
#   AGT


stopneighborsR = [unique(vcat(frameXcpairs.(stop,-x)...), dims = 1) for x in 1:2];

# Codon pairs inside an ORF, that have a stop codon in frame -2
stopNbrWithin = [stopneighborsR[x][(stopneighborsR[x][:,1] .∉ Ref(stop)) .& (stopneighborsR[x][:,2] .∉ Ref(stop)),:] for x in 1:2];

# Stop codon in frame -2, with a stop codon as an upstream codon in frame 0 #
stopNbrUPstop = [stopneighborsR[x][(stopneighborsR[x][:,1] .∈ Ref(stop)),:] for x in 1:2];

# Stop codon in frame -2, with a stop codon as an downstream codon in frame 0 #
stopNbrDNstop = [stopneighborsR[x][(stopneighborsR[x][:,1] .∉ Ref(stop)) .& (stopneighborsR2[:,2] .∈ Ref(stop)),:] for x in 1:2];

# Codons that encode a stop codon in the reverse frame
stopneighborsR0 = reverse.(comp.(stop));

p_stopWithin = [sum(nucprob.(vjoin(stopNbrWithin[x]),gccontent)) for x in 1:2];
p_stopUPstop = [sum(nucprob.(vjoin(stopNbrUPstop[x]),gccontent)) for x in 1:2];
p_stopDNstop = [sum(nucprob.(vjoin(stopNbrDNstop[x]),gccontent)) for x in 1:2];

probstopR0 = sum(nucprob.(stopneighborsR0,gccontent));

nostopneighbors = [unique(vcat(frameXcpairs.(nostopcodons,-x)...), dims = 1) for x in 1:2];

nostopNbrWithin = [nostopneighbors[x][(nostopneighbors[x][:,1] .∉ Ref(stop)) .& (nostopneighbors[x][:,2] .∉ Ref(stop)),:] for x in 1:2];

g_stopWithin = [featuregain(vjoin(nostopNbrWithin[x]),vjoin(stopNbrWithin[x]), gccontent) for x in 1:2]

l_stopWithin = [featuregain(vjoin(stopNbrWithin[x]),vjoin(nostopNbrWithin[x]), gccontent) for x in 1:2]

# Effect of selection #

function tprobsel(set1::Matrix{String},set2::Matrix{String},sdict::Dict,loss::Bool)
    s = 0
    for x in eachrow(set1)
        feasible = phcp(sdict[x[1]],sdict[x[2]])
        fset2 = set2[eachrow(set2) .∈ Ref(eachrow(feasible)),:]
        if(loss)
            s = s + featuregain([join(x)],vjoin(fset2),gccontent)
        else
            s = s + featuregain(vjoin(fset2),[join(x)],gccontent)
        end
    end
    return s
end

# Strong purifying selection #

# Stop loss #

l_stop_PurSelMax = [tprobsel(stopNbrWithin[x],nostopNbrWithin[x],syncodons,true)/p_stopWithin[x] for x in 1:2]

# Stop gain #

g_stop_PurSelMax = [tprobsel(stopNbrWithin[x],nostopNbrWithin[x],syncodons,false)/sum(nucprob.(vjoin(nostopNbrWithin[x]),gccontent)) for x in 1:2]

# Relaxed purifying selection #

# Stop loss #

l_stop_PurSelRel = [tprobsel(stopNbrWithin[x],nostopNbrWithin[x],simcodons,true)/p_stopWithin[x] for x in 1:2]

# Stop gain #

g_stop_PurSelRel = [tprobsel(stopNbrWithin[x],nostopNbrWithin[x],simcodons,false)/sum(nucprob.(vjoin(nostopNbrWithin[x]),gccontent)) for x in 1:2]
# De novo sequence divergence #

allaacodonpairs = phcp(nostopcodons,nostopcodons);

csetR = [allaacodonpairs[eachrow(allaacodonpairs) .∉ Ref(eachrow(stopNbrWithin[x])),:] for x in 1:2];

cset0 = nostopcodons[nostopcodons .∉ Ref(stopneighborsR0)]

Rcodon = (codonpair,x) -> reverse(comp(join(codonpair)[3-x+1:6-x]))

cdn2aanum = (cdn) -> aanames[transnucs(cdn)]

function aadivframeR(frm::Int,sdict::Dict)
    d = 0
    p = 0
    for Z in eachrow(csetR[frm])
        stZ = phcp(vcat(Z[1],sdict[Z[1]]),vcat(Z[2],sdict[Z[2]]))
        rcdn = Rcodon(Z,frm)
        stZ2 = stZ[(eachrow(stZ) .!= Ref(Z)) .& ([Rcodon(y,frm) for y in eachrow(stZ)] .∉ Ref(stop)) ,:]
        if(isempty(stZ2))
            dZ=0
        else
            dZ = sum([PBMEC2[cdn2aanum(rcdn),cdn2aanum(Rcodon(y,frm))] for y in eachrow(stZ2)])
        end
        pZ = nucprob(join(Z),gccontent)
        d += dZ*pZ
        p += pZ
    end
    if p==0
        return 0
    else
        return d/p
    end
end

divergence_PurSel_Max = [aadivframeR(x,syncodons) for x in 1:2]
divergence_PurSel_Rel = [aadivframeR(x,simcodons) for x in 1:2]


# # INTRONS #

# # Intron frame frequency (dmel-6.46) #

# intfrm2 = 0.423237
# intfrm1 = 0.281180
# intfrm0 = 0.295583

# # Intron nucleotide frequency (dmel-6.46) #
# first6 = [
#     0.000763	0.998862	0.000264	0.000083;
#     0.000403	0.000222	0.990159	0.009147;
#     0.601574	0.332121	0.052466	0.013755;
#     0.740860	0.097563	0.095633	0.065819;
#     0.087000	0.821628	0.065236	0.025997;
#     0.133678	0.066222	0.681662	0.118299;
# ];

# last6 = [
#     0.000347	0.998612	0.000555	0.000444;
#     0.998820	0.000194	0.000472	0.000458;
#     0.062405	0.002637	0.252797	0.682106;
#     0.268023	0.281431	0.288274	0.162188;
#     0.126433	0.047942	0.621853	0.203689;
#     0.119868	0.063307	0.588166	0.228576;
# ];

