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

nfmat = zeros(4,1);

nfmat[nucnames['A']] = 0.256362;
nfmat[nucnames['T']] = 0.204217;
nfmat[nucnames['G']] = 0.267809;
nfmat[nucnames['A']] = 0.271612;

function nucprobss(s,nfmat)
    return prod([nfmat[nucnames[x]] for x in s])
end

allcodons = kmers(3)
nostopcodons = nostop(allcodons);

noATG = setdiff(allcodons, ["ATG"]);
function ATGprobs(gccontent)
    ATGprob = nucprob("ATG",gccontent)
    ATGgain = featuregain(noATG,["ATG"],gccontent)
    ATGloss = featuregain(["ATG"],noATG,gccontent)/ATGprob
    ATGstay = featurestay(["ATG"],gccontent)
    return [ATGprob, ATGgain, ATGloss, ATGstay]
end

## Stop codon

function stopprobs(gccontent)
    stopprob = sum([nucprob(x,gccontent) for x in stopvars])
    stopgain = featuregain(nostopcodons,stopvars,gccontent)
    stoploss = featuregain(stopvars,nostop,gccontent)/stopprob
    stopstay = featurestay(stopvars,gccontent)
    return [stopprob, stopgain, stoploss, stopstay]
end

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



# p_synsubs = [featuregain([x],syncodons[x],gccontent) for x in keys(syncodons)]
# p_simsubs = [featuregain([x],simcodons[x],gccontent) for x in keys(simcodons)]

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

# Codon pairs inside an ORF, that have a stop codon in frame -x
stopNbrWithin = [stopneighborsR[x][(stopneighborsR[x][:,1] .∉ Ref(stop)) .& (stopneighborsR[x][:,2] .∉ Ref(stop)),:] for x in 1:2];

# Codons that encode a stop codon in the reverse frame
stopneighborsR0 = reverse.(comp.(stop));

# Codon pairs inside an ORF, that do not have a stop codon in frame -x
nostopneighbors = [unique(vcat(frameXcpairs.(nostopcodons,-x)...), dims = 1) for x in 1:2];

# Codons that do not encode a stop codon in the reverse frame
nostopNbr0 = nostopcodons[nostopcodons .∉ Ref(stopneighborsR0)];

nostopNbrWithin = [nostopneighbors[x][(nostopneighbors[x][:,1] .∉ Ref(stop)) .& (nostopneighbors[x][:,2] .∉ Ref(stop)),:] for x in 1:2];


function stopprobsoverlap(gccontent)
    # Probability of finding a stop codon

    p_stopWithin = [sum(nucprob.(vjoin(stopNbrWithin[x]),gccontent)) for x in 1:2];

    p_stopR0 = sum(nucprob.(stopneighborsR0,gccontent));
    push!(p_stopWithin,p_stopR0);

    # Probability of gaining a stop codon

    g_stopWithin = [featuregain(vjoin(nostopNbrWithin[x]),vjoin(stopNbrWithin[x]), gccontent) for x in 1:2]

    push!(g_stopWithin,featuregain(nostopNbr0,stopneighborsR0,gccontent));

    # Probability of losing a stop codon (conditional)

    l_stopWithin = [featuregain(vjoin(stopNbrWithin[x]),vjoin(nostopNbrWithin[x]), gccontent) for x in 1:2]

    push!(l_stopWithin,featuregain(stopneighborsR0,nostopNbr0,gccontent));

    l_stopWithin = l_stopWithin./p_stopWithin;

    s_stopWithin = [featurestay(vjoin(stopNbrWithin[x]), gccontent) for x in 1:2];

    push!(s_stopWithin,featurestay(stopneighborsR0,gccontent));


    return(hcat(p_stopWithin,g_stopWithin,l_stopWithin, s_stopWithin))
end




# Effect of selection #

function tprobsel(set1::Matrix{String},set2::Matrix{String},sdict::Dict,gccontent)
    g = 0
    l = 0
    for x in eachrow(set1)
        feasible = phcp(vcat(x[1],sdict[x[1]]),vcat(x[2],sdict[x[2]]))
        fset2 = set2[(eachrow(set2) .∈ Ref(eachrow(feasible))) .& (eachrow(set2) .!= Ref(x)),:]
        g = g + featuregain([join(x)],vjoin(fset2),gccontent)
        l = l + featuregain(vjoin(fset2),[join(x)],gccontent)
    end
    return hcat(g,l)
end

function staysel(set1::Matrix{String},sdict::Dict,gccontent)
    s = 0
    lenseq = 6
    for x in eachrow(set1)
        xj = join(x)
        feasible = phcp(vcat(x[1],sdict[x[1]]),vcat(x[2],sdict[x[2]]))
        fset1 = feasible[eachrow(feasible) .∈ Ref(eachrow(set1)),:];
        if(isempty(fset1))
            mp = 0
        else
            mp = sum([mprob(xj,y) for y in vjoin(fset1)])
        end
        s = s + nucprob(xj,gccontent)*(mp + (1-µA(mutrate,nsub))^numAT(xj)*(1-µG(mutrate,nsub))^(lenseq-numAT(xj)))
    end
    return s
end

function tprobsel0(set1::Vector{String},set2::Vector{String},sdict::Dict,gccontent)
    g = 0 
    l = 0
    for x in set1
        feasible = sdict[x]
        fset2 = set2[set2 .∈ Ref(feasible)]
        g = g + featuregain([x],fset2,gccontent)
        l = l + featuregain(fset2,[x],gccontent)
    end
    return hcat(g,l)
end

function staysel0(set1::Vector{String},sdict::Dict,gccontent)
    s = 0
    for x in set1
        feasible = sdict[x]
        fset1 = feasible[feasible .∈ Ref(set1)]
        
        if(isempty(fset1))
            mp = 0
        else
            mp = sum([mprob(x,y) for y in fset1])
        end

        s = s + nucprob(x,gccontent)*(mp + (1-µA(mutrate,nsub))^numAT(x)*(1-µG(mutrate,nsub))^(lenseq-numAT(x)))
    end
    return s
end

# Strong purifying selection #

# Stop loss #

function pstopsel(sdict,stopvals,gccontent)
    pstops = stopvals[:,1];

    GLF = vcat([tprobsel(stopNbrWithin[x],nostopNbrWithin[x],sdict,gccontent) for x in 1:2]...)

    GL0 = tprobsel0(stopneighborsR0,nostopNbr0,sdict,gccontent)

    gain = vcat(GLF[:,1],GL0[1])
    loss = vcat(GLF[:,2],GL0[2])./pstops
    stay = zeros(3,1)
    for i = 1:2
        stay[i] = staysel(stopNbrWithin[i],sdict,gccontent)
    end
    stay[3] = staysel0(stopneighborsR0,sdict,gccontent)

    return hcat(pstops,gain,loss,stay)
end

# No selection

stopvalsOL = stopprobsoverlap(gccontent)

# Relaxed purifying selection #

stopvalsOL_RelPurSel = pstopsel(simcodons,stopvalsOL,gccontent)

# Stringent purifying selection 

stopvalsOL_MaxPurSel = pstopsel(syncodons,stopvalsOL,gccontent)



function orfprobs(ATG,stop,k,ptype::Bool,polya)

    if(!ptype)
        stopprobc = stop[xProb] - polya[xProb]
        stopgainc = stop[xGain]*(1-stopprobc)/(1-stop[xProb])
    else
        stopgainc = stop[xGain]
        stopprobc = stop[xProb]
    end

    nostopstay = 1 - stopprobc - stopgainc

    orfprob = ATG[xProb]*stop[xProb]*(1 - stopprobc)^(k-2)

    gainmech1 = stop[xStay]*ATG[xGain]*nostopstay^(k-2);
    gainmech2 = ATG[xStay]*stop[xGain]*nostopstay^(k-2);
    gainmech3 = (k-2)*stop[xStay]*ATG[xStay]*stopprobc*stop[xLoss]*nostopstay^(k-1)

    orfgain =  gainmech1 + gainmech2 + gainmech3

    orfloss = stop[xLoss] + ATG[xLoss] + (k-2)*stopgainc/(1-stopprobc)

    orfstay = ATG[xStay]*stop[xStay]*(nostopstay)^(k-2)
    return [orfprob; orfgain; orfloss; orfstay]
end


# De novo sequence divergence #

allaacodonpairs = phcp(nostopcodons,nostopcodons);

csetR = [allaacodonpairs[eachrow(allaacodonpairs) .∉ Ref(eachrow(stopNbrWithin[x])),:] for x in 1:2];

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

function aadivframe0(sdict::Dict)
    d = 0
    p = 0
    for Z in nostopNbr0
        stZ = sdict[Z]
        stZ2 = stZ[stZ .!= Z]
        if(isempty(stZ2))
            dZ=0
        else
            dZ = sum([PBMEC2[cdn2aanum(Z),cdn2aanum(y)] for y in stZ2])
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
push!(divergence_PurSel_Max,aadivframe0(syncodons))

divergence_PurSel_Rel = [aadivframeR(x,simcodons) for x in 1:2]
push!(divergence_PurSel_Rel,aadivframe0(simcodons))


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

