using Plots
using Measures
cd(Base.source_dir())
include("nucleotidefuncts.jl")

organism = "dmel"

cm2pt = (cm) -> 28.3465*cm
figdir = joinpath(Base.source_dir(),"../Manuscripts/Figures/M2_main/");
colors = ["#FFCC00","#5599FF","#D40000","#754473","#000000"];
lstyles = [:solid,:dash,:dot]

default(linecolor = :black, linewidth = 2, tickfont = font(10,"Helvetica"), 
guidefont = font(13,"Helvetica"),framestyle = :box, legend = false);

stopcodons = ["TAA","TAG","TGA"]
function nostop(codonset)
	codonset[codonset .∉ Ref(stopcodons)]
end

(xProb, xGain, xLoss, xStay) = [i for i in 1:4];

allcodons = kmers(3)
nostopcodons = nostop(allcodons);

nsub, nucsmbt = readdlm(organism*"_mutbias.txt",'\t',header = true);
if organism == "scer"
    mutrate = 1.7e-10
end

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
    stopprob = sum([nucprob(x,gccontent) for x in stopcodons])
    stopgain = featuregain(nostopcodons,stopcodons,gccontent)
    stoploss = featuregain(stopcodons,nostopcodons,gccontent)/stopprob
    stopstay = featurestay(stopcodons,gccontent)
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

allaacodonpairs = phcp(nostopcodons,nostopcodons);

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


stopneighborsR = [unique(vcat(frameXcpairs.(stopcodons,-x)...), dims = 1) for x in 1:2];

# Codon pairs inside an ORF, that have a stop codon in frame -x
stopNbrWithin = [stopneighborsR[x][(stopneighborsR[x][:,1] .∉ Ref(stopcodons)) .& (stopneighborsR[x][:,2] .∉ Ref(stopcodons)),:] for x in 1:2];

# Codons that encode a stop codon in the reverse frame
stopneighborsR0 = reverse.(comp.(stopcodons));

# Codon pairs inside an ORF, that do not have a stop codon in frame -x
nostopneighbors = [unique(vcat(frameXcpairs.(nostopcodons,-x)...), dims = 1) for x in 1:2];

# Codons that do not encode a stop codon in the reverse frame
nostopNbr0 = nostopcodons[nostopcodons .∉ Ref(stopneighborsR0)];

nostopNbrWithin = [nostopneighbors[x][(nostopneighbors[x][:,1] .∉ Ref(stopcodons)) .& (nostopneighbors[x][:,2] .∉ Ref(stopcodons)),:] for x in 1:2];

csetR = [allaacodonpairs[eachrow(allaacodonpairs) .∉ Ref(eachrow(stopNbrWithin[x])),:] for x in 1:2];

Rcodon = (codonpair,x) -> reverse(comp(join(codonpair)[3-x+1:6-x]))
cdn2aanum = (cdn) -> aanames[transnucs(cdn)]


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
    lenseq = 3
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

function orfprobs(ATG,stop,k)

    stopgainc = stop[xGain]
    stopprobc = stop[xProb]

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

ncodons = [30:300;];
gcrange = [0.3:0.1:0.6;];

orfvalsITG = zeros(length(gcrange),length(ncodons),4);
orfvalsONS, orfvalsORS, orfvalsOSS = [zeros(length(gcrange),length(ncodons),3,4) for i = 1:3];

PstpR = zeros(3,length(gcrange));
GstpR, LstpR = [zeros(3,3,length(gcrange)) for i = 1:2];

for g in eachindex(gcrange)
    atgvals = ATGprobs(gcrange[g])

    # Intergenic region
    stopitg = stopprobs(gcrange[g])

    # No selection
    stopvalsOL = stopprobsoverlap(gcrange[g])
    PstpR[:,g] = stopvalsOL[:,1]/stopitg[1]

    # Relaxed purifying selection #
    stopvalsOL_RelPurSel = pstopsel(simcodons,stopvalsOL,gcrange[g])
    
    # Stringent purifying selection 
    stopvalsOL_MaxPurSel = pstopsel(syncodons,stopvalsOL,gcrange[g])

    GstpR[:,:,g] = hcat(stopvalsOL[:,2], stopvalsOL_RelPurSel[:,2], stopvalsOL_RelPurSel[:,2])/stopitg[2];

    LstpR[:,:,g] = hcat(stopvalsOL[:,3], stopvalsOL_RelPurSel[:,3], stopvalsOL_RelPurSel[:,3])/stopitg[3];

    for k in eachindex(ncodons)
        orfvalsITG[g,k,:] = orfprobs(atgvals,stopitg,ncodons[k])
        for f in 1:3
            orfvalsONS[g,k,f,:] = orfprobs(atgvals,stopvalsOL[f,:],ncodons[k])
            orfvalsORS[g,k,f,:] = orfprobs(atgvals,stopvalsOL_RelPurSel[f,:],ncodons[k])
            orfvalsOSS[g,k,f,:] = orfprobs(atgvals,stopvalsOL_MaxPurSel[f,:],ncodons[k])
        end
    end
    println("Done: ",gcrange[g])
end

include("antisenseGenes_supplement.jl")

orfvalsITG_X = zeros(length(ncodons),4);
orfvalsONS_X, orfvalsORS_X, orfvalsOSS_X = [zeros(length(ncodons),3,4) for i = 1:3];

atgvals_X = ATGprobsX(trimerfreq)
stopitg_X = stopprobsX(trimerfreq)
stopvalsONS_X = stopprobsoverlapX(codonfreq,dicodonfreq)
stopvalsORS_X = pstopselX(simcodons,stopvalsONS_X,codonfreq,dicodonfreq)
stopvalsOSS_X = pstopselX(syncodons,stopvalsONS_X,codonfreq,dicodonfreq)

for k in eachindex(ncodons)
    orfvalsITG_X[k,:] = orfprobs(atgvals_X,stopitg_X,ncodons[k])
    for f in 1:3
        orfvalsONS_X[k,f,:] = orfprobs(atgvals_X,stopvalsONS_X[f,:],ncodons[k])
        orfvalsORS_X[k,f,:] = orfprobs(atgvals_X,stopvalsORS_X[f,:],ncodons[k])
        orfvalsOSS_X[k,f,:] = orfprobs(atgvals_X,stopvalsOSS_X[f,:],ncodons[k])
    end
end

pORF = plot(xlabel = "ORF length (codons)", ylabel = "Relative ORF Probability\n (Frame 1)", size = (width = cm2pt(11), height = cm2pt(10)));

for q in 1:4
    plot!(pORF,ncodons,log2.(orfvalsONS[q,:,1,1]./orfvalsITG[q,:,1]),
        linecolor = colors[q]
    );
end

savefig(pORF, figdir*"pORF_antisense_"*organism*"GC.pdf")

pORFX = plot(xlabel = "ORF length (codons)", ylabel = "Relative ORF Probability\n log2(antisense/intergenic)", size = (width = cm2pt(11), height = cm2pt(10)));

for q = 1:3
    plot!(pORFX,ncodons,log2.(orfvalsONS_X[:,q,1]./orfvalsITG_X[:,1]),
            linestyle = lstyles[q]
    );
end

savefig(pORFX, figdir*"pORF_antisense_"*organism*"WC.pdf")

plots_gain, plots_loss = [Array{Plots.Plot{Plots.GRBackend}}(undef,3,3) for i = 1:2]; 

rat = zeros(length(ncodons),4,3,2);
gx = [1:4;];

pnames = ["stationary","gain","loss"];
snames = ["No selection", "Relaxed", "Stringent"];

for f = 1:3
    for p = 1:2
        rat[:,1:4,1,p] = log2.(orfvalsONS[gx,:,f,p+1]./orfvalsITG[gx,:,p+1])'
        rat[:,1:4,2,p] = log2.(orfvalsORS[gx,:,f,p+1]./orfvalsITG[gx,:,p+1])'
        rat[:,1:4,3,p] = log2.(orfvalsOSS[gx,:,f,p+1]./orfvalsITG[gx,:,p+1])'
    end
    for s in 1:3
        ydatag = rat[:,:,s,1]
        lbg = minimum(ydatag); ubg = maximum(ydatag)
        mng = maximum(abs.(ydatag));
        mng>10 ? d = 2 : d = 3
        plots_gain[f,s] = plot(ncodons, ydatag[:,1],
            linecolor = colors[1],
            yticks = round.(range(lbg,stop=ubg,length = 4),digits=2)
        );

        ydatal = rat[:,:,s,2]
        lbl = minimum(ydatal); ubl = maximum(ydatal)
        mnl = maximum(abs.(ydatal));
        mnl>10 ? d = 2 : d = 3
        plots_loss[f,s] = plot(ncodons, ydatal[:,1],
            linecolor = colors[1],
            yticks = round.(range(lbl,stop=ubl,length = 4),digits=2)
        );
        for q = 2:4
            plot!(plots_gain[f,s],ncodons,ydatag[:,q], linecolor = colors[q]);
            plot!(plots_loss[f,s],ncodons,ydatal[:,q], linecolor = colors[q]);
        end

    end
end

pG = plot(plots_gain..., size = (width = cm2pt(27.5), height = cm2pt(25)));
savefig(pG, figdir*"pORFgain_antisense_"*organism*"GC.pdf")

pL = plot(plots_loss..., size = (width = cm2pt(27.5), height = cm2pt(25)));
savefig(pL, figdir*"pORFloss_antisense_"*organism*"GC.pdf")

plots_gain_dmel, plots_loss_dmel = [Array{Plots.Plot{Plots.GRBackend}}(undef,1,3) for i = 1:2];


rat_dmel = zeros(length(ncodons),3,2);

for f = 1:3
    for p = 1:2
        rat_dmel[:,1,p] = log2.(orfvalsONS_X[:,f,p+1]./orfvalsITG_X[:,p+1])'
        rat_dmel[:,2,p] = log2.(orfvalsORS_X[:,f,p+1]./orfvalsITG_X[:,p+1])'
        rat_dmel[:,3,p] = log2.(orfvalsOSS_X[:,f,p+1]./orfvalsITG_X[:,p+1])'
    end

    ydatag = rat_dmel[:,:,1]
    lbg = minimum(ydatag); ubg = maximum(ydatag)
    mng = maximum(abs.(ydatag));
    mng>10 ? d = 2 : d = 3
    plots_gain_dmel[f] = plot(ncodons, ydatag[:,1],
        linestyle = lstyles[1],
        yticks = round.(range(lbg,stop=ubg,length = 4),digits=1)
    );

    ydatal = rat_dmel[:,:,2]
    lbl = minimum(ydatal); ubl = maximum(ydatal)
    mnl = maximum(abs.(ydatal));
    mnl>10 ? d = 2 : d = 3
    plots_loss_dmel[f] = plot(ncodons, ydatal[:,1],
        linestyle = lstyles[1],
         yticks = round.(range(lbl,stop=ubl,length = 4),digits=1)
    );

    for q = 2:3
        plot!(plots_gain_dmel[f],ncodons,ydatag[:,q], linestyle = lstyles[q]);
        plot!(plots_loss_dmel[f],ncodons,ydatal[:,q], linestyle = lstyles[q]);
    end

end

pGD = plot(plots_gain_dmel..., layout = (1,3), size = (width = cm2pt(27.5), height = cm2pt(9)));
savefig(pGD, figdir*"pORFgain_antisense_"*organism*"WC.pdf")

pLD = plot(plots_loss_dmel..., layout = (1,3), size = (width = cm2pt(27.5), height = cm2pt(9)));
savefig(pLD, figdir*"pORFloss_antisense_"*organism*"WC.pdf")

# De novo sequence divergence #


function aadivframeR(frm::Int,sdict::Dict,gccontent)
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
        pZ = nucprob(join(Z),gccontent)
        d += abs(dZ)*pZ
        p += pZ
    end
    if p==0
        return 0
    else
        return d/p
    end
end

function aadivframe0(sdict::Dict,gccontent)
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
        pZ = nucprob(Z,gccontent)
        d += abs(dZ)*pZ
        p += pZ
    end
    if p==0
        return 0
    else
        return d/p
    end
end

divergence_PurSel_Max, divergence_PurSel_Rel = [zeros(length(gcrange),3) for i = 1:2];

for g in eachindex(gcrange)
    gcx = gcrange[g]
    divergence_PurSel_Max[g,1:2] = [aadivframeR(x,syncodons,gcx) for x in 1:2]
    divergence_PurSel_Max[g,3] = aadivframe0(syncodons,gcx)

    divergence_PurSel_Rel[g,1:2] = [aadivframeR(x,simcodons,gcx) for x in 1:2]
    divergence_PurSel_Rel[g,3] = aadivframe0(simcodons,gcx)
end

divergence_PurSel_Rel_X = [aadivframeR_X(x,simcodons,dicodonfreq) for x in 1:2]
push!(divergence_PurSel_Rel_X,aadivframe0_X(simcodons,codonfreq))

divergence_PurSel_Max_X = [aadivframeR_X(x,syncodons,dicodonfreq) for x in 1:2]
push!(divergence_PurSel_Max_X,aadivframe0_X(syncodons,codonfreq))

alldivRel = vcat(divergence_PurSel_Rel,divergence_PurSel_Rel_X')
alldivMax = vcat(divergence_PurSel_Max,divergence_PurSel_Max_X')

bwf = 1/7;

pDivRel = plot(
    xlabel = "Frame",
    ylabel = "Divergence", 
    size = (width = cm2pt(11.5), height = cm2pt(10)),   
);


for j = 1:3
    bar!(pDivRel,[j+x*bwf for x in 1:5],alldivRel[:,j],
        fill = colors,
        bar_width = bwf,
        linecolor = nothing
    )
end

xticks!(pDivRel,[1:3;] .+ 3*bwf, string.([1,2,0]));

savefig(pDivRel,figdir*"DivergenceRel_"*organism*".pdf")


pDivMax = plot(
    xlabel = "Frame",
    ylabel = "Divergence", 
    size = (width = cm2pt(11.5), height = cm2pt(10)),   
);

for j = 1:3
    bar!(pDivMax,[j+x*bwf for x in 1:5],alldivMax[:,j],
        fill = colors,
        bar_width = bwf,
        linecolor = nothing
    )
end

xticks!(pDivMax,[1:3;] .+ 3*bwf, string.([1,2,0]));

savefig(pDivMax,figdir*"DivergenceMax_"*organism*".pdf")