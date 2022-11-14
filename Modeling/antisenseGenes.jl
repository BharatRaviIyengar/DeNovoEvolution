cd(Base.source_dir())
include("nucleotidefuncts.jl")
stop = ["TAA","TAG","TGA"]

function nostop(codonset)
	codonset[codonset .∉ Ref(stop)]
end

gccontent = 0.53942 # all CDS #

p_T = 0.204217;
p_A = 0.256362;
p_G = 0.267809;
p_C = 0.271612;

allcodons = kmers(3)
nostopcodons = nostop(allcodons);

codonpairprob = (cpair) -> sum(nucprob.(cpair[:,1],gccontent) .* nucprob.(cpair[:,2],gccontent))

syncodons = Dict(x => gencode[(gencode[:,3].== transnucs(x)) .& (gencode[:,1] .!= x),1] for x in nostopcodons)

PBMEC = 0 .+readdlm("PMBEC.txt")[2:21,2:21]
aanames = Dict(aas[i] => i for i in 1:20)

simcodons = syncodons
for x in keys(simcodons)
    y = transnucs(x)
    simaas = aas[PBMEC[:,aanames[y]] .> 0.05]
    simaas = simaas[simaas .!= y]
    simcody = gencode[gencode[:,3] .∈ Ref(simaas),1]
    simcodons[x] = vcat(simcodons[x],simcody)
end

p_synsubs = [featuregain([x],syncodons[x],gccontent) for x in keys(syncodons)]
p_simsubs = [featuregain([x],simcodons[x],gccontent) for x in keys(simcodons)]

# Codon pairs that have a stop codon in frame -2
# Example: 
# TTCACG
#  AGT
stopneighborsR2 = unique(vcat(frameXcpairs.(stop,-2)...), dims = 1);

# Codon pairs inside an ORF, that have a stop codon in frame -2
stopNbrR2within = stopneighborsR2[(stopneighborsR2[:,1] .∉ Ref(stop)) .& (stopneighborsR2[:,2] .∉ Ref(stop)),:];

# Stop codon in frame -2, with a stop codon as an upstream codon in frame 0 #
stopNbrR2stopUP = stopneighborsR2[(stopneighborsR2[:,1] .∈ Ref(stop)),:];

# Stop codon in frame -2, with a stop codon as an downstream codon in frame 0 #
stopNbrR2stopDN = stopneighborsR2[(stopneighborsR2[:,1] .∉ Ref(stop)) .& (stopneighborsR2[:,2] .∈ Ref(stop)),:];

# Stop codon in frame -2, with a start codon as an upstream codon in frame 0 #
stopNbrR2initUP = stopneighborsR2[(stopneighborsR2[:,1] .== "ATG"),:];

# Stop codon in frame -2, with a start codon as an downstream codon in frame 0 #
stopNbrR2initDN = stopneighborsR2[(stopneighborsR2[:,1] .== "ATG") .& (stopneighborsR2[:,2].∈ Ref(stop)),:];

# Codon pairs that have a stop codon in frame -1
# Example: 
# TCTCAG
#   AGT
stopneighborsR1 = unique(vcat(frameXcpairs.(stop,-1)...), dims = 1);

# Stop codon in frame -1, with a stop codon as an upstream codon in frame 0 #
stopNbrR1stopUP = stopneighborsR1[(stopneighborsR1[:,1] .∈ Ref(stop)),:];

# Stop codon in frame -1, with a stop codon as an downstream codon in frame 0 #
stopNbrR1stopDN = stopneighborsR1[(stopneighborsR1[:,1] .∉ Ref(stop)) .& (stopneighborsR1[:,2] .∈ Ref(stop)),:];

# Codon pairs inside an ORF, that have a stop codon in frame -1
stopNbrR1within = stopneighborsR1[(stopneighborsR1[:,1] .∉ Ref(stop)) .& (stopneighborsR1[:,2] .∉ Ref(stop)),:];

# Stop codon in frame -1, with a start codon as an upstream codon in frame 0 #
stopNbrR2initUP = stopneighborsR2[(stopneighborsR1[:,1] .== "ATG"),:];

# Stop codon in frame -1, with a start codon as an downstream codon in frame 0 #
stopNbrR2initDN = stopneighborsR2[(stopneighborsR1[:,2] .== "ATG") .& (stopneighborsR2[:,2].∈ Ref(stop)),:];

# Codons that encode a stop codon in the reverse frame
stopneighborsR0 = reverse.(comp.(stop));

p_stopR2within = codonpairprob(stopNbrR2within);
p_stopR2stopUP = codonpairprob(stopNbrR2stopUP);
p_stopR2stopDN = codonpairprob(stopNbrR2stopDN);
p_stopR2initUP = codonpairprob(stopNbrR2initUP);
p_stopR2initDN = codonpairprob(stopNbrR2initDN);

p_stopR1within = codonpairprob(stopNbrR1within);
p_stopR1stopUP = codonpairprob(stopNbrR1stopUP);
p_stopR1stopDN = codonpairprob(stopNbrR1stopDN);
p_stopR1initUP = codonpairprob(stopNbrR1initUP);
p_stopR1initDN = codonpairprob(stopNbrR1initDN);

probstopR0 = sum(nucprob.(stopneighborsR0,gccontent));

initneighborsR1 = unique(frameXcpairs("ATG",-1), dims = 1);
initNbrR1within = initneighborsR1[(initneighborsR1[:,1] .∉ Ref(stop)) .& (initneighborsR1[:,2].∉ Ref(stop)),:];

initNbrR1stopUP = initneighborsR1[(initneighborsR1[:,1] .∈ Ref(stop)),:];

initNbrR1stopDN = initneighborsR1[(initneighborsR1[:,1] .∉ Ref(stop)) .& (initneighborsR1[:,2] .∈ Ref(stop)),:];

initneighborsR2 = unique(frameXcpairs("ATG",-2), dims = 1);
initNbrR2within = initneighborsR2[(initneighborsR2[:,1] .∉ Ref(stop)) .& (initneighborsR2[:,2].∉ Ref(stop)),:];

initNbrR2stopUP = initneighborsR2[(initneighborsR2[:,1] .∈ Ref(stop)),:];

initNbrR2stopDN = initneighborsR2[(initneighborsR2[:,1] .∉ Ref(stop)) .& (initneighborsR2[:,2] .∈ Ref(stop)),:];

probinitR1 = codonpairprob(initNbrR1within);
probinitR2 = codonpairprob(initNbrR1within);


# INTRONS #

# Intron frame frequency (dmel-6.46) #

intfrm2 = 0.423237
intfrm1 = 0.281180
intfrm0 = 0.295583

# Intron nucleotide frequency (dmel-6.46) #
first6 = [
    0.000763	0.998862	0.000264	0.000083;
    0.000403	0.000222	0.990159	0.009147;
    0.601574	0.332121	0.052466	0.013755;
    0.740860	0.097563	0.095633	0.065819;
    0.087000	0.821628	0.065236	0.025997;
    0.133678	0.066222	0.681662	0.118299;
];

last6 = [
    0.000347	0.998612	0.000555	0.000444;
    0.998820	0.000194	0.000472	0.000458;
    0.062405	0.002637	0.252797	0.682106;
    0.268023	0.281431	0.288274	0.162188;
    0.126433	0.047942	0.621853	0.203689;
    0.119868	0.063307	0.588166	0.228576;
];

