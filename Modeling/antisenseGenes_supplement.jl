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

nprob3 = x -> codonfreq[codonfreq[:,1] .== x,2]
nprob6 = x -> dicodonfreq[dicodonfreq[:,1] .== x,2]

function featuregain3(set1,set2)
    if(isempty(set1) || isempty(set2))
        @warn "empty feature sets, returning 0"
        return 0
    end
    v1 = sum([nprob3(s1)*sum([mprob(s1,s2) for s2 in set2])
                for s1 in set1])
    return v1
end

function featurestay3(set1)
    v =  sum(
        [nprob3(s1)*(
            sum([mprob(s1,s2) for s2 in set1]) +
            (1-µA(mutrate,nsub))^numAT(s1)*
            (1-µG(mutrate,nsub))^(3-numAT(s1))
            ) for s1 in set1]
        )
    return v
end

function featuregain6(set1,set2)
    if(isempty(set1) || isempty(set2))
        @warn "empty feature sets, returning 0"
        return 0
    end
    v1 = sum([nprob6(s1)*sum([mprob(s1,s2) for s2 in set2])
                for s1 in set1])
    return v1
end

function featurestay6(set1)
    v =  sum(
        [nprob6(s1)*(
            sum([mprob(s1,s2) for s2 in set1]) +
            (1-µA(mutrate,nsub))^numAT(s1)*
            (1-µG(mutrate,nsub))^(6-numAT(s1))
            ) for s1 in set1]
        )
    return v
end

noATG = setdiff(allcodons, ["ATG"]);
function ATGprobsX(gccontent)
    ATGprob = nprob3("ATG",gccontent)
    ATGgain = featuregain3(noATG,["ATG"])
    ATGloss = featuregain3(["ATG"],noATG)/ATGprob
    ATGstay = featurestay3(["ATG"])
    return [ATGprob, ATGgain, ATGloss, ATGstay]
end

## Stop codon

function stopprobsX(gccontent)
    stopprob = sum([nucprob(x,gccontent) for x in stopcodons])
    stopgain = featuregain3(nostopcodons,stopcodons)
    stoploss = featuregain3(stopcodons,nostopcodons)/stopprob
    stopstay = featurestay3(stopcodons)
    return [stopprob, stopgain, stoploss, stopstay]
end

function stopprobsoverlapX()
    # Probability of finding a stop codon

    p_stopWithin = [sum(nprob6.(vjoin(stopNbrWithin[x]),gccontent)) for x in 1:2];

    p_stopR0 = sum(nprob6.(stopneighborsR0));
    push!(p_stopWithin,p_stopR0);

    # Probability of gaining a stop codon

    g_stopWithin = [featuregain6(vjoin(nostopNbrWithin[x]),vjoin(stopNbrWithin[x])) for x in 1:2]

    push!(g_stopWithin,featuregain6(nostopNbr0,stopneighborsR0));

    # Probability of losing a stop codon (conditional)

    l_stopWithin = [featuregain6(vjoin(stopNbrWithin[x]),vjoin(nostopNbrWithin[x])) for x in 1:2]

    push!(l_stopWithin,featuregain6(stopneighborsR0,nostopNbr0));

    l_stopWithin = l_stopWithin./p_stopWithin;

    s_stopWithin = [featurestay6(vjoin(stopNbrWithin[x])) for x in 1:2];

    push!(s_stopWithin,featurestay(stopneighborsR0,gccontent));


    return(hcat(p_stopWithin,g_stopWithin,l_stopWithin, s_stopWithin))
end




# Effect of selection #

function tprobselX(set1::Matrix{String},set2::Matrix{String},sdict::Dict)
    g = 0
    l = 0
    for x in eachrow(set1)
        feasible = phcp(vcat(x[1],sdict[x[1]]),vcat(x[2],sdict[x[2]]))
        fset2 = set2[(eachrow(set2) .∈ Ref(eachrow(feasible))) .& (eachrow(set2) .!= Ref(x)),:]
        g = g + featuregain6([join(x)],vjoin(fset2))
        l = l + featuregain6(vjoin(fset2),[join(x)])
    end
    return hcat(g,l)
end

function stayselX(set1::Matrix{String},sdict::Dict)
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
        s = s + nprob6(xj)*(mp + (1-µA(mutrate,nsub))^numAT(xj)*(1-µG(mutrate,nsub))^(lenseq-numAT(xj)))
    end
    return s
end

function tprobsel0X(set1::Vector{String},set2::Vector{String},sdict::Dict)
    g = 0 
    l = 0
    for x in set1
        feasible = sdict[x]
        fset2 = set2[set2 .∈ Ref(feasible)]
        g = g + featuregain6([x],fset2)
        l = l + featuregain6(fset2,[x])
    end
    return hcat(g,l)
end

function staysel0X(set1::Vector{String},sdict::Dict)
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

        s = s + nprob6(x)*(mp + (1-µA(mutrate,nsub))^numAT(x)*(1-µG(mutrate,nsub))^(lenseq-numAT(x)))
    end
    return s
end

# Strong purifying selection #

# Stop loss #

function pstopselX(sdict,stopvals)
    pstops = stopvals[:,1];

    GLF = vcat([tprobselX(stopNbrWithin[x],nostopNbrWithin[x],sdict) for x in 1:2]...)

    GL0 = tprobsel0X(stopneighborsR0,nostopNbr0,sdict)

    gain = vcat(GLF[:,1],GL0[1])
    loss = vcat(GLF[:,2],GL0[2])./pstops
    stay = zeros(3,1)
    for i = 1:2
        stay[i] = stayselX(stopNbrWithin[i],sdict)
    end
    stay[3] = staysel0X(stopneighborsR0,sdict)

    return hcat(pstops,gain,loss,stay)
end