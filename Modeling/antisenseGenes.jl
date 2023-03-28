using Plots
using Measures
cd(Base.source_dir())
include("nucleotidefuncts.jl")

organism = "scer"

cm2pt = (cm) -> 28.3465*cm
figdir = joinpath(Base.source_dir(),"../Manuscripts/Figures/M2_main/pdf/");
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
    p_stopR12 = [sum(nucprob.(vjoin(stopNbrWithin[x]),gccontent)) for x in 1:2];
    p_stopR0 = sum(nucprob.(stopneighborsR0,gccontent));
    p_stopWithin = vcat(p_stopR0,p_stopR12);

    # Probability of gaining a stop codon

    g_stopR12 = [featuregain(vjoin(nostopNbrWithin[x]),vjoin(stopNbrWithin[x]), gccontent) for x in 1:2]
    g_stopR0 = featuregain(nostopNbr0,stopneighborsR0,gccontent);
    g_stopWithin = vcat(g_stopR0,g_stopR12);

    # Probability of losing a stop codon (conditional)

    l_stopR12 = [featuregain(vjoin(stopNbrWithin[x]),vjoin(nostopNbrWithin[x]), gccontent) for x in 1:2]
    l_stopR0 = featuregain(stopneighborsR0,nostopNbr0,gccontent);
	l_stopWithin = vcat(l_stopR0,l_stopR12)./p_stopWithin; 

    s_stopR12 = [featurestay(vjoin(stopNbrWithin[x]), gccontent) for x in 1:2];
    s_stopR0 = featurestay(stopneighborsR0,gccontent);
    s_stopWithin = vcat(s_stopR0,s_stopR12);

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

    gain = vcat(GL0[1],GLF[:,1])
    loss = vcat(GL0[2],GLF[:,2])./pstops
    stay = zeros(3,1)
    stay[1] = staysel0(stopneighborsR0,sdict,gccontent)
    for i = 2:3
        stay[i] = staysel(stopNbrWithin[i-1],sdict,gccontent)
    end
    

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
    ubg = maximum(rat[:,:,:,1]);
    mng = maximum(abs.(rat[:,:,:,1]));
    lbg = minimum(rat[:,:,:,1]);
    mng>10 ? dg = 2 : dg = 3
    yrngg = round.(range(lbg,stop=ubg,length = 4),digits=2);
    ylg = [lbg - 0.1*(ubg - lbg), ubg + 0.1*(ubg -lbg)]

    ubl = maximum(rat[:,:,:,2]);
    mnl = maximum(abs.(rat[:,:,:,2]));
    lbl = minimum(rat[:,:,:,2]);
    mnl>10 ? dl = 2 : dl = 3
    yrngl = round.(range(lbl,stop=ubl,length = 4),digits=2);
    yll = [lbl - 0.05*(ubl -lbl), ubl + 0.05*(ubl -lbl)]
    for s in 1:3
        ydatag = rat[:,:,s,1]
        # lbg = minimum(ydatag); ubg = maximum(ydatag)
        # mng = maximum(abs.(ydatag));
        # mng>10 ? d = 2 : d = 3
        plots_gain[s,f] = plot(ncodons, ydatag[:,1],
            linecolor = colors[1],
            ylims = ylg,
            yticks = yrngg,
            xticks = [80, 160, 240]
        );

        ydatal = rat[:,:,s,2]
        # lbl = minimum(ydatal); ubl = maximum(ydatal)
        # mnl = maximum(abs.(ydatal));
        # mnl>10 ? d = 2 : d = 3
        plots_loss[s,f] = plot(ncodons, ydatal[:,1],
            linecolor = colors[1],
            ylims = yll,
            yticks = yrngl,
            xticks = [80, 160, 240]
        );
        for q = 2:4
            plot!(plots_gain[s,f],ncodons,ydatag[:,q], linecolor = colors[q]);
            plot!(plots_loss[s,f],ncodons,ydatal[:,q], linecolor = colors[q]);
        end

    end
end

pG = plot(plots_gain..., size = (width = cm2pt(21), height = cm2pt(18)));
savefig(pG, figdir*"pORFgain_antisense_"*organism*"GC.pdf")

pL = plot(plots_loss..., size = (width = cm2pt(21), height = cm2pt(18)));
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
        yticks = round.(range(lbg,stop=ubg,length = 4),digits=1),
        xticks = [80, 160, 240]
    );

    ydatal = rat_dmel[:,:,2]
    lbl = minimum(ydatal); ubl = maximum(ydatal)
    mnl = maximum(abs.(ydatal));
    mnl>10 ? d = 2 : d = 3
    plots_loss_dmel[f] = plot(ncodons, ydatal[:,1],
        linestyle = lstyles[1],
         yticks = round.(range(lbl,stop=ubl,length = 4),digits=1),
         xticks = [80, 160, 240]
    );

    for q = 2:3
        plot!(plots_gain_dmel[f],ncodons,ydatag[:,q], linestyle = lstyles[q]);
        plot!(plots_loss_dmel[f],ncodons,ydatal[:,q], linestyle = lstyles[q]);
    end

end

pGD = plot(plots_gain_dmel..., layout = (1,3), size = (width = cm2pt(21), height = cm2pt(6.75)));
savefig(pGD, figdir*"pORFgain_antisense_"*organism*"WC.pdf")

pLD = plot(plots_loss_dmel..., layout = (1,3), size = (width = cm2pt(21), height = cm2pt(6.75)));
savefig(pLD, figdir*"pORFloss_antisense_"*organism*"WC.pdf")

# De novo sequence divergence #


function aadivframeR(frm::Int,sdict::Dict,gccontent)
    pvec,dvec = [zeros(size(csetR[frm],1),1) for i = 1:2]
    ind = 0
    for Z in eachrow(csetR[frm])
        ind += 1
            stZ = phcp(vcat(Z[1],sdict[Z[1]]),vcat(Z[2],sdict[Z[2]]))
        rcdn = Rcodon(Z,frm)
        stZ2 = stZ[(eachrow(stZ) .!= Ref(Z)) .& ([Rcodon(y,frm) for y in eachrow(stZ)] .∉ Ref(stopcodons)) ,:]
        if(isempty(stZ2))
            dZ = 0
            mZ = 0
        else
            dZZ = [PBMEC2[cdn2aanum(rcdn),cdn2aanum(Rcodon(y,frm))] for y in eachrow(stZ2)]

            mZZ = [mprob(rcdn,Rcodon(y,frm)) for y in eachrow(stZ2)]

            dZ = sum(dZZ.*mZZ)
            mZ = sum(mZZ)
        end
        pZ = nucprob(join(Z),gccontent)
        dvec[ind] = abs(dZ)*pZ
        pvec[ind] = pZ*mZ
    end
    p = sum(pvec)
    dm = sum(dvec)/p
    dv = sum((dvec - dm.*pvec).^2)/p
    if p==0
        return [0,0]
    else
        return [dm,sqrt(dv)]
    end
end

function aadivframe0(sdict::Dict,gccontent)
    pvec,dvec = [zeros(size(nostopNbr0)) for i = 1:2]
    ind = 0
    for Z in nostopNbr0
        ind += 1
        stZ = sdict[Z]
        stZ2 = stZ[stZ .!= Z]
        if(isempty(stZ2))
            dZ = 0
            mZ = 0
        else
            dZZ = [PBMEC2[cdn2aanum(Z),cdn2aanum(y)] for y in stZ2]
            mZZ = [mprob(Z,y) for y in stZ2]

            dZ = sum(dZZ.*mZZ)
            mZ = sum(mZZ)
        end
        pZ = nucprob(Z,gccontent)
        dvec[ind] = abs(dZ)*pZ
        pvec[ind] = pZ*mZ
    end
    p = sum(pvec)
    dm = sum(dvec)/p
    dv = sum((dvec - dm.*pvec).^2)/p
    if p==0
        return [0,0]
    else
        return [dm,sqrt(dv)]
    end
end

function divall(gcc)
    pvec,dvec = [zeros(size(nostopcodons)) for i = 1:2]
    ind = 0
    for Z in nostopcodons
        ind += 1
        dZZ = [PBMEC2[cdn2aanum(Z),cdn2aanum(y)] for y in nostopcodons]
        mZZ = [mprob(Z,y) for y in nostopcodons]

        dZ = sum(dZZ.*mZZ)
        mZ = sum(mZZ)

        pZ = nucprob(Z,gcc)
        dvec[ind] = abs(dZ)*pZ
        pvec[ind] = pZ*mZ
    end
    p = sum(pvec)
    dm = sum(dvec)/p
    dv = sum((dvec - dm.*pvec).^2)/p
    if p==0
        return [0,0]
    else
        return [dm,sqrt(dv)]
    end
end

divergence_PurSel_Max, divergence_PurSel_Rel = [zeros(length(gcrange),3,2) for i = 1:2];

divergence_PurSel_Max_X, divergence_PurSel_Rel_X = [zeros(3,2) for i = 1:2];

for g in eachindex(gcrange)
    gcx = gcrange[g]
    for x = 1:2
        divergence_PurSel_Max[g,x+1,:] = aadivframeR(x,syncodons,gcx)
        divergence_PurSel_Rel[g,x+1,:] = aadivframeR(x,simcodons,gcx)
    end
    divergence_PurSel_Max[g,1,:] = aadivframe0(syncodons,gcx)
    divergence_PurSel_Rel[g,1,:] = aadivframe0(simcodons,gcx)
end

divergence_PurSel_Rel_X[1,:] = aadivframe0_X(simcodons,codonfreq);
divergence_PurSel_Max_X[1,:] = aadivframe0_X(syncodons,codonfreq);
for x = 1:2
    divergence_PurSel_Rel_X[x+1,:] = aadivframeR_X(x,simcodons,dicodonfreq)
    divergence_PurSel_Max_X[x+1,:] = aadivframeR_X(x,syncodons,dicodonfreq)
end


alldivRel_avg = vcat(divergence_PurSel_Rel[:,:,1],divergence_PurSel_Rel_X[:,1]')
alldivMax_avg = vcat(divergence_PurSel_Max[:,:,1],divergence_PurSel_Max_X[:,1]')

alldivRel_std = vcat(divergence_PurSel_Rel[:,:,2],divergence_PurSel_Rel_X[:,2]')
alldivMax_std = vcat(divergence_PurSel_Max[:,:,2],divergence_PurSel_Max_X[:,2]')

alldivNoS = hcat([divall(x) for x in gcrange]...)

bwf = 1/7;

pDivRel = plot(
    xlabel = "Frame",
);
pDivMax = plot(
    xlabel = "Frame",
);

for j = 1:3
    bar!(pDivRel,[j+x*bwf for x in 1:5],alldivRel_avg[:,j],
        fill = colors,
        bar_width = bwf,
        linecolor = nothing
    )
    plot!(pDivRel,[j+x*bwf for x in 1:5],alldivNoS[1,:],
        markersize = 3,
        markershape = :circle,
        markercolor = colors,
        markerstrokecolor = :white,
        markerstrokewidth = 2.5,
        linecolor = nothing
    )
    bar!(pDivMax,[j+x*bwf for x in 1:5],alldivMax_avg[:,j],
        fill = colors,
        bar_width = bwf,
        linecolor = nothing
    )
    plot!(pDivMax,[j+x*bwf for x in 1:5],alldivNoS[1,:],
        markersize = 3,
        markershape = :circle,
        markercolor = colors,
        markerstrokecolor = :white,
        markerstrokewidth = 2.5,
        linecolor = nothing
    )
end

xticks!(pDivMax,[1:3;] .+ 3*bwf, string.([0,1,2]));
xticks!(pDivRel,[1:3;] .+ 3*bwf, string.([0,1,2]));

ylims!(pDivRel,ylims(pDivMax));
yticks!(pDivRel,yticks(pDivMax)[1][1]);

pDivAll = plot(pDivRel,pDivMax, layout = (1,2), size = (width = cm2pt(20), height = cm2pt(10)))

savefig(pDivAll,figdir*"Divergence_"*organism*".pdf")