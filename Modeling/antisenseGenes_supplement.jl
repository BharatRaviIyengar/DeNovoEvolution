using DelimitedFiles

normalize = x -> x/sum(x)

codonfreq = readdlm(joinpath(Base.source_dir(),"codons.txt"), '\t');
codonfreq[codonfreq[:,1] .∈ Ref(stopcodons),2] .= 0;
codonfreq[:,2] = codonfreq[:,2]/sum(codonfreq[:,2]);

dcfreqdata = readdlm(joinpath(Base.source_dir(),"dicodons.txt"), '\t');
cc = z -> [x[z:z+2] for x in dcfreqdata[:,1]]

dcfreqdata[(cc(1) .∈ Ref(stopcodons)) .| (cc(3) .∈ Ref(stopcodons)) ,2] .= 0;

dicodonfreq = hcat(kmers(6),zeros(4^6,1));
for x in eachrow(dcfreqdata)
    dicodonfreq[dicodonfreq[:,1] .== x[1],2] .= x[2];
end

dicodonfreq[:,2] = normalize(dicodonfreq[:,2]);

trimerfreq = readdlm(joinpath(Base.source_dir(),"trimers2.txt"), '\t');
trimerfreq[:,2] = normalize(trimerfreq[:,2])

function ATGprobsX(t3)
    ATGprob = nprob3("ATG",t3)
    ATGgain = featuregain3(noATG,["ATG"],t3)
    ATGloss = featuregain3(["ATG"],noATG,t3)/ATGprob
    ATGstay = featurestay3(["ATG"],t3)
    return [ATGprob, ATGgain, ATGloss, ATGstay]
end

## Stop codon

function stopprobsX(t3)
    stopprob = sum([nprob3(x,t3) for x in stopcodons])
    stopgain = featuregain3(nostopcodons,stopcodons,t3)
    stoploss = featuregain3(stopcodons,nostopcodons,t3)/stopprob
    stopstay = featurestay3(stopcodons,t3)
    return [stopprob, stopgain, stoploss, stopstay]
end

function stopprobsoverlapX(t3,t6)
    # Probability of finding a stop codon

    p_stopWithin = [sum([nprob6(y,t6) for y in vjoin(stopNbrWithin[x])]) for x in 1:2];

    p_stopR0 = sum([nprob3(y,t3) for y in stopneighborsR0]);
    push!(p_stopWithin,p_stopR0);

    # Probability of gaining a stop codon

    g_stopWithin = [featuregain6(vjoin(nostopNbrWithin[x]),vjoin(stopNbrWithin[x]),t6) for x in 1:2]

    push!(g_stopWithin,featuregain3(nostopNbr0,stopneighborsR0,t3));

    # Probability of losing a stop codon (conditional)

    l_stopWithin = [featuregain6(vjoin(stopNbrWithin[x]),vjoin(nostopNbrWithin[x]),t6) for x in 1:2]

    push!(l_stopWithin,featuregain3(stopneighborsR0,nostopNbr0,t3));

    l_stopWithin = l_stopWithin./p_stopWithin;

    s_stopWithin = [featurestay6(vjoin(stopNbrWithin[x]),t6) for x in 1:2];

    push!(s_stopWithin,featurestay3(stopneighborsR0,t3));

    return(hcat(p_stopWithin,g_stopWithin,l_stopWithin, s_stopWithin))
end




# Effect of selection #

function tprobselX(set1::Matrix{String},set2::Matrix{String},sdict::Dict,t6)
    g = 0
    l = 0
    for x in eachrow(set1)
        feasible = phcp(vcat(x[1],sdict[x[1]]),vcat(x[2],sdict[x[2]]))
        fset2 = set2[(eachrow(set2) .∈ Ref(eachrow(feasible))) .& (eachrow(set2) .!= Ref(x)),:]
        g = g + featuregain6([join(x)],vjoin(fset2),t6)
        l = l + featuregain6(vjoin(fset2),[join(x)],t6)
    end
    return hcat(g,l)
end

function stayselX(set1::Matrix{String},sdict::Dict,t6)
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
        s = s + nprob6(xj,t6)*(mp + (1-µA(mutrate,nsub))^numAT(xj)*(1-µG(mutrate,nsub))^(lenseq-numAT(xj)))
    end
    return s
end

function tprobsel0X(set1::Vector{String},set2::Vector{String},sdict::Dict,t3)
    g = 0 
    l = 0
    for x in set1
        feasible = sdict[x]
        fset2 = set2[set2 .∈ Ref(feasible)]
        g = g + featuregain3([x],fset2,t3)
        l = l + featuregain3(fset2,[x],t3)
    end
    return hcat(g,l)
end

function staysel0X(set1::Vector{String},sdict::Dict,t3)
    s = 0
    lenseq = 3
    for x in set1
        feasible = sdict[x]
        fset1 = feasible[feasible .∈ Ref(set1)]
        
        if(isempty(fset1))
            mp = 0
        else
            mp = sum([mprob(x,y) for y in fset1])
        end

        s = s + nprob3(x,t3)*(mp + (1-µA(mutrate,nsub))^numAT(x)*(1-µG(mutrate,nsub))^(lenseq-numAT(x)))
    end
    return s
end

# Strong purifying selection #

# Stop loss #

function pstopselX(sdict,stopvals,t3,t6)
    pstops = stopvals[:,1];

    GLF = vcat([tprobselX(stopNbrWithin[x],nostopNbrWithin[x],sdict,t6) for x in 1:2]...)

    GL0 = tprobsel0X(stopneighborsR0,nostopNbr0,sdict,t3)

    gain = vcat(GLF[:,1],GL0[1])
    loss = vcat(GLF[:,2],GL0[2])./pstops
    stay = zeros(3,1)
    for i = 1:2
        stay[i] = stayselX(stopNbrWithin[i],sdict,t6)
    end
    stay[3] = staysel0X(stopneighborsR0,sdict,t3)

    return hcat(pstops,gain,loss,stay)
end

function aadivframeR_X(frm::Int,sdict::Dict,t6)
    d = 0
    p = 0
    for Z in eachrow(csetR[frm])
        stZ = phcp(vcat(Z[1],sdict[Z[1]]),vcat(Z[2],sdict[Z[2]]))
        rcdn = Rcodon(Z,frm)
        stZ2 = stZ[(eachrow(stZ) .!= Ref(Z)) .& ([Rcodon(y,frm) for y in eachrow(stZ)] .∉ Ref(stopcodons)) ,:]
        if(isempty(stZ2))
            dZ=0
        else
            dZ = sum([PBMEC2[cdn2aanum(rcdn),cdn2aanum(Rcodon(y,frm))] for y in eachrow(stZ2)])
        end
        pZ = nprob6(join(Z),t6)
        d += abs(dZ)*pZ
        p += pZ
    end
    if p==0
        return 0
    else
        return d/p
    end
end

function aadivframe0_X(sdict::Dict,t3)
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
        pZ = nprob3(join(Z),t3)
        d += abs(dZ)*pZ
        p += pZ
    end
    if p==0
        return 0
    else
        return d/p
    end
end