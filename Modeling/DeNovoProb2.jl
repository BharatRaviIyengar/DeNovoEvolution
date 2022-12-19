using Combinatorics
using Plots
using Measures
# using HypothesisTests

include("nucleotidefuncts.jl")

cm2pt = (cm) -> 28.3465*cm
figdir = joinpath(Base.source_dir(),"../Manuscripts/Figures/");

default(linecolor = :black, linewidth = 2, tickfont = font(10,"Helvetica"), 
guidefont = font(13,"Helvetica"),framestyle = :box, legend = false);

aaorder = ["W","F","Y","I","V","L","M","C","D","E","G","A","P","H","K","R","S","T","N","Q"];
ixaaorder = indexin(aaorder,aas);
gccontent = 0.42

(xProb, xGain, xLoss, xStay) = [i for i in 1:4];

mers6 = kmers(6);
allcodons = kmers(3)

function tataprobs(gccontent)
    ## TATA box
    # Consensus sequence = TATAWAWR
    # is TSS - 25 - 8nt (transcript must be +33nt longer)
    tataboxes = vec(join.(product(("TATA" .* collect("TA") .* "A"), "TA", "GA")))

    tatavars = unique(reduce(vcat,
        reduce(vcat, [
            [degen(s,1,nucs) for s in tataboxes],
            [degen(s,2,"TA") for s in tataboxes],
            [degen(s,3,"TAC") for s in tataboxes],
            [degen(s,4,"AT") for s in tataboxes],
            [degen(s,6,nucs) for s in tataboxes],
            [degen(s,7,nucs) for s in tataboxes],
            [degen(s,8,nucs) for s in tataboxes],
        ])
    ));

    notata = setdiff(kmers(8),tatavars);

    tataprob = sum([nucprob(x,gccontent) for x in tatavars]);
    tatagain = featuregain(notata,tatavars,gccontent);
    tataloss = featuregain(tatavars,notata,gccontent)/tataprob;
    tatastay = featurestay(tatavars,gccontent);
    return [tataprob,tatagain,tataloss,tatastay]
end

function inrprobs(gccontent)
    ## Initiator element
    # Consensus sequence = "BBCABW";  DOI:10.1101/gad.295980.117
    inrvars = vec(join.(product("TGC","TGC","C", "A","TGC","TA")));
    noinrs = setdiff(mers6, inrvars);
    inrprob = sum([nucprob(x,gccontent) for x in inrvars]);
    inrgain = featuregain(noinrs,inrvars,gccontent);
    inrloss = featuregain(inrvars,noinrs,gccontent);
    inrstay = featurestay(inrvars,gccontent)/inrprob;
    return [inrprob,inrgain,inrloss,inrstay]
end

function polyaprobs(gccontent)
    ## PolyA signal
    polyavars = ["AATAAA";"ATTAAA"; "AGTAAA";"TATAAA"];
    nopolyas = setdiff(mers6, polyavars);
    polyaprob = sum([nucprob(x,gccontent) for x in polyavars]);
    polyagain = featuregain(nopolyas,polyavars,gccontent);
    polyaloss = featuregain(polyavars,nopolyas,gccontent)/polyaprob;
    polyastay = featurestay(polyavars,gccontent);
    return [polyaprob,polyagain,polyaloss,polyastay]
end

function rnaprobs(tata,inr,polya,orflen,ptype::Bool)
    
    # Number of possible polyA signals given ORF is present
    # TAA within polyA signal cannot exist with intact ORF
    if(ptype)
        noPAsites = 3*orflen - 5 # total 6-mers
    else
        noPAsites = 2*orflen - 3 # 6-mers not in with ORF
    end
    nopolyaprob = 1 - polya[xProb]
    nopolyastay = nopolyaprob - polya[xGain]

    ## Transcription
    # (Initiator or TATA) and polyA

    # Probability of finding a transcript of length ≥ orflen
    rnaprob = (tata[xProb] + inr[xProb]) * polya[xProb] * nopolyaprob^noPAsites;
    
    # RNA gain mechanism 1
    # Gain of any promoter
    # polyA stays at the end of sequence
    # no polyA present in the sequence 
    gainmech1 = (tata[xGain] + inr[xGain]) * polya[xStay]*nopolyastay^noPAsites 

    # RNA gain mechanism 2
    # Any one promoter stays
    # polyA gained at the end of sequence
    # no polyA present in the sequence 
    gainmech2 = (tata[xStay] + inr[xStay]) * polya[xGain] * nopolyastay^noPAsites

    # RNA gain mechanism 3
    # Any one promoter stays
    # polyA stays at the end of sequence
    # A polyA present in the sequence is lost (for every 6mer in the sequence)
    gainmech3 = (noPAsites-1)*(tata[xStay] + inr[xStay]) * polya[xStay] * polya[xLoss] * polya[xProb] * nopolyastay^(noPAsites-1);


    rnagain =  gainmech1 + gainmech2 + gainmech3
    
    # Fraction of Inr-only promoters
    Qinr = inr[xProb]*(1-tata[xProb])/(tata[xProb] + inr[xProb] - inr[xProb]*tata[xProb]);

    # Fraction of TATA-only promoters
    Qtata = tata[xProb]*(1-inr[xProb])/(tata[xProb] + inr[xProb] - inr[xProb]*tata[xProb]);

    # Fraction of promoters containing both TATA and Inr
    Qboth = 1 - Qtata - Qinr;
    notatastay  = 1 - tata[xGain] - tata[xProb];
    noinrstay = 1 - inr[xGain] - inr[xProb]
    rnaloss = (Qinr*inr[xLoss]*notatastay +
      Qtata*tata[xLoss]*noinrstay +
      Qboth*tata[xLoss]*inr[xLoss] +
      polya[xLoss] + noPAsites*polya[xGain]
    );

    rnastay = (inr[xStay] + tata[xStay])*polya[xStay]*nopolyastay^noPAsites;
    return [rnaprob; rnagain; rnaloss; rnastay]
end

function ATGprobs(gccontent)
    ## Start codon
    noATG = setdiff(allcodons, ["ATG"])
    ATGprob = nucprob("ATG",gccontent)
    ATGgain = featuregain(noATG,["ATG"],gccontent)
    ATGloss = featuregain(["ATG"],noATG,gccontent)/ATGprob
    ATGstay = featurestay(["ATG"],gccontent)
    return [ATGprob, ATGgain, ATGloss, ATGstay]
end

function stopprobs(gccontent)
    ## Stop codon
    stopvars = ["TAA","TAG","TGA"]
    nostop = setdiff(allcodons, stopvars)
    stopprob = sum([nucprob(x,gccontent) for x in stopvars])
    stopgain = featuregain(nostop,stopvars,gccontent)
    stoploss = featuregain(stopvars,nostop,gccontent)/stopprob
    stopstay = featurestay(stopvars,gccontent)
    return [stopprob, stopgain, stoploss, stopstay]
end



function orfprobs(ATG,stop,k,ptype::Bool,polya)

    if(!ptype)
        stop[xProb] -= polya[xProb]
        stop[xStay] -= polya[x]
    end
    nostopstay = 1 - stop[xProb] - stop[xGain]
    orfprob = ATG[xProb]*stop[xProb]*(1 - stop[xProb])^(k-2)

    orfgain = nostopstay^(k-2)*(
    stop[xStay]*ATG[xGain] + ATG[xStay]*stop[xGain] +
    (k-2)*(stop[xStay]*ATG[xStay]*stop[xProb]*stop[xLoss]/nostopstay)
    )
    orfloss = stop[xLoss] + ATG[xLoss] + (k-2)*stop[xGain]/(1-stop[xProb])
    orfstay = ATG[xStay]*stop[xStay]*(nostopstay)^(k-2)
    return [orfprob; orfgain; orfloss; orfstay]
end

ncodons = [30:300;];

ATGvals = ATGprobs(gccontent);
stopvals = stopprobs(gccontent);
tatavals = tataprobs(gccontent);
inrvals = inrprobs(gccontent);
polyavals = polyaprobs(gccontent);

(rnaprob, rnagain, rnaloss, rnastay) = [zeros(size(ncodons)) for i = 1:4 ];
(crnaprob, crnagain, crnaloss, crnastay) = [zeros(size(ncodons)) for i = 1:4 ];
(orfprob, orfgain, orfloss, orfstay) = [zeros(size(ncodons)) for i = 1:4 ];

for k in eachindex(ncodons)
    (orfprob[k], orfgain[k], orfloss[k], orfstay[k]) = orfprobs(ATGvals,stopvals,ncodons[k])
    (rnaprob[k], rnagain[k], rnaloss[k], rnastay[k]) = rnaprobs(tatavals, inrvals, polyavals, ncodons[k],true)
    (crnaprob[k], crnagain[k], crnaloss[k], crnastay[k]) = rnaprobs(tatavals, inrvals, polyavals, ncodons[k],false)
end

genegain = ((rnagain .+ rnastay) .* orfgain + orfstay .* crnagain);
genegain2 = genegain./(1 .-rnaprob .- orfprob);
geneloss = orfloss .+ rnaloss;

rnafirst = rnastay .* orfgain;
orffirst = orfstay .* crnagain;

onlyrnagain = (1 .- orfprob .- orfgain).*rnagain;
onlyorfgain = (1 .- rnaprob .- rnagain).*orfgain;

rnafirst2 = onlyrnagain.*(1 .- rnaloss).*(orfgain./(1 .-orfprob));
orffirst2 = onlyorfgain.*(1 .- orfloss).*(rnagain./(1 .-rnaprob));

plot_genegain_geneloss = plot(ncodons,log.(genegain2)./log.(geneloss),
    xlabel = "ORF length (codons)",
    ylabel = "# Gene Losses \n per Gene Gain",
    yticks = [2:0.4:3.2;],
    size = (width = cm2pt(12), height = cm2pt(11)),
);
savefig(plot_genegain_geneloss, figdir*"geneGainLoss_new.pdf")

gcrange = [0.3:0.01:0.6;];

(rnaprobGC, rnagainGC, rnalossGC, rnastayGC) = [zeros(length(gcrange),length(ncodons)) for i = 1:4 ];
(crnaprobGC, crnagainGC, crnalossGC, crnastayGC) = [zeros(length(gcrange),length(ncodons)) for i = 1:4 ];
(orfprobGC, orfgainGC, orflossGC, orfstayGC) = [zeros(length(gcrange),length(ncodons)) for i = 1:4 ];

for g in eachindex(gcrange)
    ATGvalsG = ATGprobs(gcrange[g]);
    stopvalsG = stopprobs(gcrange[g]);
    tatavalsG = tataprobs(gcrange[g]);
    inrvalsG = inrprobs(gcrange[g]);
    polyavalsG = polyaprobs(gcrange[g]);
    for k in eachindex(ncodons)
        (orfprobGC[g,k], orfgainGC[g,k], orflossGC[g,k], orfstayGC[g,k]) = orfprobs(ATGvalsG,stopvalsG,ncodons[k])
        (rnaprobGC[g,k], rnagainGC[g,k], rnalossGC[g,k], rnastayGC[g,k]) = rnaprobs(tatavalsG, inrvalsG, polyavalsG, ncodons[k],true)
        (crnaprobGC[g,k], crnagainGC[g,k], crnalossGC[g,k], crnastayGC[g,k]) = rnaprobs(tatavalsG, inrvalsG, polyavalsG, ncodons[k],false)
    end
    println("Done: ",g)
end

k50 = findfirst(ncodons .==200)

genegain50 = ((rnagainGC[:,k50] .+ rnastayGC[:,k50]) .* orfgainGC[:,k50] .+ orfstayGC[:,k50] .* crnagainGC[:,k50]);
genegain2_50 = genegain50./(1 .-rnaprobGC[:,k50] .- orfprobGC[:,k50]);
geneloss50 = orflossGC[:,k50] .+ rnalossGC[:,k50]

# genegainX = ((rnagainGC .+ rnastayGC) .* orfgainGC .+ orfstayGC .* crnagainGC);
# genegainX2 = genegainX./(1 .-rnaprobGC .- orfprobGC);
# genelossX = orflossGC .+ rnalossGC

# genegainlossratio = log.(genegainX2)./log.(genelossX)

# whichGCloss = [findfirst(x.>=2) for x in eachcol(genegainlossratio)];

plot_genegain_geneloss_50_gcrange = plot(gcrange,log.(genegain2_50)./log.(geneloss50),
    xlabel = "GC-content",
    ylabel = "# Gene Losses \n per Gene Gain",
    size = (width = cm2pt(12.4), height = cm2pt(11)),
);
savefig(plot_genegain_geneloss_50_gcrange, figdir*"geneGainLoss50gcr_new.pdf")

plot_rnafist_orffirst = plot(ncodons,log2.(orffirst./rnafirst),
    xlabel = "ORF length (codons)",
    ylabel = "P_{ORF-first}\nP_RNA-first",
    size = (width = cm2pt(12.5), height = cm2pt(11)),
);
savefig(plot_rnafist_orffirst, figdir*"whoisfirst_new.pdf")

rnafirstgcr = hcat([rnastayGC[x,:] .* orfgainGC[x,:] for x in eachindex(gcrange)]...);
orffirstgcr =  hcat([orfstayGC[x,:] .* crnagainGC[x,:] for x in eachindex(gcrange)]...);

whichfirstgcr = [findfirst(x.<=0) for x in eachcol(log2.(orffirstgcr./rnafirstgcr))];

whichfirstgcr[isnothing.(whichfirstgcr)] .= 1;

plot_whichfirstgcr = plot(gcrange,ncodons[whichfirstgcr],
    xlabel = "GC-content",
    ylabel = "Minimum codons for \n RNA-first trajectory",
    size = (width = cm2pt(13), height = cm2pt(11)),
);

savefig(plot_whichfirstgcr, figdir*"whoisfirstgcr_new.pdf")

plot_rnafist_orffirst2 = plot(ncodons[1:50],log.(orffirst2./rnafirst2)[1:50],
    xlabel = "ORF length (codons)",
    ylabel = "P_{ORF-first}\nP_RNA-first",
    size = (width = cm2pt(12), height = cm2pt(11)),
    xticks = [35:7:80;]
);
savefig(plot_rnafist_orffirst2, figdir*"first_ORF_RNA2_new.pdf")

tl = log10.(rnaloss./orfgain);
ol = log10.(orfloss/rnagain);

plot_Loss = plot(ncodons,log2.(orfloss./rnaloss),
    xlabel = "ORF length (codons)",
    ylabel = "P_ORF-loss\nP_RNA-loss",
    size = (width = cm2pt(12), height = cm2pt(11)),
);

savefig(plot_Loss, figdir*"pLoss_new.pdf")

plot_loss_50_gcrange = plot(gcrange,log2.(orflossGC[:,k50]./rnalossGC[:,k50]),
    xlabel = "GC-content",
    ylabel = "P_ORF-loss\nP_RNA-loss",
    size = (width = cm2pt(12), height = cm2pt(11)),
);
savefig(plot_loss_50_gcrange, figdir*"pLoss50gcr_new.pdf")

onlyrnaloss = orfstay .*(rnaprob*rnaloss)

aaprob = [sum(nucprob.(gencode[gencode[:,3] .== x,1],gccontent)) for x in aas];
aaprob2 = aaprob./sum(aaprob);
exp_hydropathy = sum(hydropathy3.*aaprob2)

acidic = indexin(["D","E"],aas);
basic = indexin(["K","R"],aas);
hydrophobic = hydropathy3 .<0;

allcodons = kmers(3);
nostopcodons = allcodons[allcodons .∉ Ref(["TAA","TAG","TGA"])]

pHydrophobic = sum(aaprob[hydrophobic]);

plot_aafreq = bar(aaorder,aaprob2[ixaaorder],
    xticks = (0.5:20, aaorder),
    yticks = (0:0.03:0.9),
    fill = :black,
    xlabel = "Amino acid",
    ylabel = "Expected frequency\n(Normalized)",
    size = (width = cm2pt(15), height = cm2pt(11))
);
savefig(plot_aafreq, figdir*"pAAfreq.pdf")


aasubprob = zeros(20,20);

for i = 1:20
        for j = 1:20
                if(i==j)
                        continue
                end
                set1 = gencode[gencode[:,3] .== aas[i],1];
                set2 = gencode[gencode[:,3] .== aas[j],1];
                aasubprob[i,j] = featuregain(set1,set2,gccontent);
        end
end

aasubhighprob = aasubprob.*(aasubprob .>mutrate^2)

aasubprob_conditional = aasubhighprob./aaprob2;

asymmetric_transitions = log2.(aasubhighprob./aasubhighprob');
asymmetric_conditional = log2.(aasubprob_conditional./aasubprob_conditional');

# c = kmeans(asymmetric_transitions,3);
# idx = sortperm(assignments(c));
# hcl = hclust(aasubhighprob .+ aasubhighprob');
# idx = hcl.order;

gainAA = Vector{Float64}(vec(sum(aasubprob, dims=1)'));
lossAA = Vector{Float64}(vec(sum(aasubprob, dims=2)./aaprob2));
stayAA = [sum(featurestay(gencode[gencode[:,3] .== x,1],gccontent)) for x in aas];

getHP = (a) -> hydropathy3[aas .== a]

codonshydrophobic = gencode[(in).(gencode[:,3],Ref(aas[hydrophobic])),1];

codonshydrophilic = kmers(3)[(!in).(kmers(3),Ref(codonshydrophobic))];

hmp = 0
mp=0
for i in nostopcodons
    for j in nostopcodons
        x = mprob(i,j)
        yj = getHP(transnucs(j))[1]
        yi = getHP(transnucs(i))[1]
        if(sign(yi)!=sign(yj))
            y = yj - yi
        else
            y = 0
        end
        z = nucprob(i,gccontent)
        hmp += x*y*z
        mp += x
    end
end
MutΔHydrophobic = hmp/mp


toHydrophobic = featuregain(codonshydrophilic,codonshydrophobic, gccontent);

frHydrophobic = featuregain(codonshydrophobic,codonshydrophilic, gccontent);

stayHydrophobic = featurestay(codonshydrophobic, gccontent);

plot_aasubprob = heatmap(aaorder,aaorder,log2.(aasubhighprob[ixaaorder,ixaaorder]),
    xticks = (0.5:20, aaorder), 
    yticks = (0.5:20, aaorder),
    size = (width = cm2pt(15), height = cm2pt(14)),
    colorbar_title = "Substitution Probability",
    legend = :bottom
);

savefig(plot_aasubprob, figdir*"pAASub.pdf")

plot_aasymprob = heatmap(aaorder,aaorder,asymmetric_transitions[ixaaorder,ixaaorder], 
    xticks = (0.5:20, aaorder), 
    yticks = (0.5:20,aaorder),
    size = (width = cm2pt(15), height = cm2pt(14)),
    colorbar_title = "Substitution likelihood",
    legend = :bottom,
    ylabel = "Original amino acid",
    xlabel = "Substituted amino acid",
    c = :bluesreds
);

savefig(plot_aasymprob, figdir*"pAAsubSym.pdf")

plot_aasymprob2 = heatmap(aaorder,aaorder,asymmetric_conditional[ixaaorder,ixaaorder], 
    xticks = (0.5:20, aaorder), 
    yticks = (0.5:20,aaorder),
    size = (width = cm2pt(15), height = cm2pt(14)),
    colorbar_title = "Substitution likelihood",
    legend = :bottom,
    ylabel = "Original amino acid",
    xlabel = "Substituted amino acid",
    c = :bluesreds
);
savefig(plot_aasymprob2, figdir*"pAAsubSymCond.svg")

plot_toaa = bar(aaorder,gainAA[ixaaorder]./sum(gainAA),
    xticks = (0.5:20, aaorder),
    yticks = (0:0.03:0.9),
    fill = :black,
    xlabel = "Amino acid",
    ylabel = "Gain Probability\n(Normalized)",
    size = (width = cm2pt(15), height = cm2pt(11))
);

savefig(plot_toaa, figdir*"pToAA.pdf")

plot_fraa = bar(aaorder,lossAA[ixaaorder]./sum(lossAA),
    xticks = (0.5:20, aaorder),
    fill = :black,
    xlabel = "Amino acid",
    ylabel = "Loss Probability\n(Normalized)",
    size = (width = cm2pt(15), height = cm2pt(11))
);

savefig(plot_fraa, figdir*"pFrAA.pdf")

