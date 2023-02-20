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
gccontent = 0.41

(xProb, xGain, xLoss, xStay) = [i for i in 1:4];

mers6 = kmers(6);

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

function tataprobs(gccontent)
    tataprob = sum([nucprob(x,gccontent) for x in tatavars]);
    tatagain = featuregain(notata,tatavars,gccontent);
    tataloss = featuregain(tatavars,notata,gccontent)/tataprob;
    tatastay = featurestay(tatavars,gccontent);
    return [tataprob,tatagain,tataloss,tatastay]
end

## Initiator element
# Consensus sequence = "BBCABW";  DOI:10.1101/gad.295980.117
inrvars = vec(join.(product("TGC","TGC","C", "A","TGC","TA")));
noinrs = setdiff(mers6, inrvars);

function inrprobs(gccontent)
    inrprob = sum([nucprob(x,gccontent) for x in inrvars]);
    inrgain = featuregain(noinrs,inrvars,gccontent);
    inrloss = featuregain(inrvars,noinrs,gccontent);
    inrstay = featurestay(inrvars,gccontent);
    return [inrprob,inrgain,inrloss,inrstay]
end

## PolyA signal
polyavars = ["AATAAA";"ATTAAA"; "AGTAAA";"TATAAA"];
nopolyas = setdiff(mers6, polyavars);

function polyaprobs(gccontent)
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
        noPAsites = 2*orflen - 5 # 6-mers not in with ORF
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

    rnaloss = polya[xLoss] + noPAsites*polya[xGain];

    rnastay = (inr[xStay] + tata[xStay])*polya[xStay]*nopolyastay^noPAsites;
    return [rnaprob; rnagain; rnaloss; rnastay]
end


function rnaprobs_nopromoter(polya,orflen,ptype::Bool)
    
    # Number of possible polyA signals given ORF is present
    # TAA within polyA signal cannot exist with intact ORF
    if(ptype)
        noPAsites = 3*orflen - 5 # total 6-mers
    else
        noPAsites = 2*orflen - 5 # 6-mers not in with ORF
    end
    nopolyaprob = 1 - polya[xProb]
    nopolyastay = nopolyaprob - polya[xGain]

    ## Transcription
    # (Initiator or TATA) and polyA

    # Probability of finding a transcript of length ≥ orflen
    rnaprob = polya[xProb] * nopolyaprob^noPAsites;
    
    # RNA gain mechanism 1
    # polyA gained at the end of sequence
    # no polyA present in the sequence 

    gainmech1 = polya[xGain] * nopolyastay^noPAsites

    # RNA gain mechanism 2
    # polyA stays at the end of sequence
    # A polyA present in the sequence is lost (for every 6mer in the sequence)

    gainmech2 = (noPAsites-1)* polya[xStay] * polya[xLoss] * polya[xProb] * nopolyastay^(noPAsites-1);

    rnagain =  gainmech1 + gainmech2 
    
    rnaloss = polya[xLoss] + noPAsites*polya[xGain];

    rnastay = polya[xStay]*nopolyastay^noPAsites;

    return [rnaprob; rnagain; rnaloss; rnastay]
end

function kimura(s,ne)
    fp = (1 - exp(-2*s))/(1-exp(-4*ne*s))
    return fp
end

allcodons = kmers(3);

## Start codon
noATG = setdiff(allcodons, ["ATG"]);
function ATGprobs(gccontent)
    ATGprob = nucprob("ATG",gccontent)
    ATGgain = featuregain(noATG,["ATG"],gccontent)
    ATGloss = featuregain(["ATG"],noATG,gccontent)/ATGprob
    ATGstay = featurestay(["ATG"],gccontent)
    return [ATGprob, ATGgain, ATGloss, ATGstay]
end

## Stop codon
stopvars = ["TAA","TAG","TGA"];
nostop = setdiff(allcodons, stopvars);
function stopprobs(gccontent)
    stopprob = sum([nucprob(x,gccontent) for x in stopvars])
    stopgain = featuregain(nostop,stopvars,gccontent)
    stoploss = featuregain(stopvars,nostop,gccontent)/stopprob
    stopstay = featurestay(stopvars,gccontent)
    return [stopprob, stopgain, stoploss, stopstay]
end



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

ncodons = [30:300;];

gcrange = [0.3:0.01:0.6;];

(rnaprob, rnagain, rnaloss, rnastay) = [zeros(length(gcrange),length(ncodons)) for i = 1:4 ];
(crnaprob, crnagain, crnaloss, crnastay) = [zeros(length(gcrange),length(ncodons)) for i = 1:4 ];
(orfprob, orfgain, orfloss, orfstay) = [zeros(length(gcrange),length(ncodons)) for i = 1:4 ];
(corfprob, corfgain, corfloss, corfstay) = [zeros(length(gcrange),length(ncodons)) for i = 1:4 ];

for g in eachindex(gcrange)
    ATGvalsG = ATGprobs(gcrange[g]);
    stopvalsG = stopprobs(gcrange[g]);
    # tatavalsG = tataprobs(gcrange[g]);
    # inrvalsG = inrprobs(gcrange[g]);
    polyavalsG = polyaprobs(gcrange[g]);
    for k in eachindex(ncodons)

        (orfprob[g,k], orfgain[g,k], orfloss[g,k], orfstay[g,k]) = orfprobs(ATGvalsG,stopvalsG,ncodons[k],true,polyavalsG);

        (corfprob[g,k], corfgain[g,k], corfloss[g,k], corfstay[g,k]) = orfprobs(ATGvalsG,stopvalsG,ncodons[k],false,polyavalsG);

        (rnaprob[g,k], rnagain[g,k], rnaloss[g,k], rnastay[g,k]) = rnaprobs_nopromoter(polyavalsG, ncodons[k],true)

        (crnaprob[g,k], crnagain[g,k], crnaloss[g,k], crnastay[g,k]) = rnaprobs_nopromoter(polyavalsG, ncodons[k],false)
    
    end
    println("Done: ",g)
end

promprob = 5331240/43164203; # CRMs in open chromatin / Total intergenic region

rnaprob = rnaprob*promprob;
rnagain = rnagain*promprob;
crnaprob = crnaprob*promprob;
crnagain = crnagain*promprob;

dGC = findall(gcrange .∈ Ref([0.34:0.08:0.5;]));
cxlen = [30:30:150;]
klen = ncodons .∈ Ref(cxlen)

# Gene gain/loss probabilities

genegain = rnastay .* corfgain + orfstay .* crnagain + orfgain.*rnagain;

genegain2 = genegain./(1 .-rnaprob .- orfprob);

geneloss = orfloss.+ rnaloss;

# Plot gene gain/loss probabilities versus ORF length

plot_genegain_geneloss = plot(ncodons,(log.(genegain2)./log.(geneloss))[dGC[1],:],
    xlabel = "ORF length (codons)",
    ylabel = "# Gene Losses \n per Gene Gain",
    size = (width = cm2pt(12), height = cm2pt(11)),
    linestyle = :dash,
    legend = :bottom,
    label = Int(100*gcrange[dGC[1]])
);

plot!(plot_genegain_geneloss,ncodons,(log.(genegain2)./log.(geneloss))[dGC[2],:],
    linestyle = :solid, 
    label = Int(100*gcrange[dGC[2]])
);

plot!(plot_genegain_geneloss,ncodons,(log.(genegain2)./log.(geneloss))[dGC[3],:],
    linestyle = :dot,
    label = Int(100*gcrange[dGC[3]])
);

savefig(plot_genegain_geneloss, figdir*"geneGainLoss_promoterless.pdf")

# Plot gene gain/loss probabilities versus GC content

plot_genegain_geneloss_klen_gcrange = plot(100*gcrange,(log.(genegain2)./log.(geneloss))[:,klen],
    xlabel = "GC %",
    ylabel = "# Gene Losses \n per Gene Gain",
    size = (width = cm2pt(12.5), height = cm2pt(11)),
    xlims = [29,65]
);
savefig(plot_genegain_geneloss_klen_gcrange, figdir*"geneGainLoss_klengcr_promoterless.pdf")

# Single step trajectory probability

rnafirst = rnastay .* corfgain;
orffirst = orfstay .* crnagain;

# Two step trajectory probability

onlyrnagain = (1 .- orfprob .- orfgain).*rnagain;
onlyorfgain = (1 .- rnaprob .- rnagain).*orfgain;

rnafirst2 = onlyrnagain.*(1 .- rnaloss).*(corfgain./(1 .-orfprob));
orffirst2 = onlyorfgain.*(1 .- orfloss).*(crnagain./(1 .-rnaprob));

# Plot single step trajectory probability versus ORF length

plot_rnafist_orffirst = plot(ncodons,log2.(rnafirst./orffirst)[dGC[1],:],
    xlabel = "ORF length (codons)",
    ylabel = "P_{ORF-first}\nP_RNA-first",
    size = (width = cm2pt(12.5), height = cm2pt(11)),
    linestyle = :dash,
    legend = :bottom,
    label = Int(100*gcrange[dGC[1]])
);
    

plot!(plot_rnafist_orffirst, ncodons,log2.(rnafirst./orffirst)[dGC[2],:],
    linestyle = :solid,
    label = Int(100*gcrange[dGC[2]])
);


plot!(plot_rnafist_orffirst, ncodons,log2.(rnafirst./orffirst)[dGC[3],:],
    linestyle = :dot,
    label = Int(100*gcrange[dGC[3]])
);

savefig(plot_rnafist_orffirst, figdir*"whoisfirst_promoterless.pdf")

# Plot two step trajectory probability versus ORF length

plot_rnafist_orffirst2 = plot(ncodons,log2.(rnafirst2./orffirst2)[dGC[1],:],
    xlabel = "ORF length (codons)",
    ylabel = "P_{ORF-first}\nP_RNA-first",
    size = (width = cm2pt(12.5), height = cm2pt(11)),
    linestyle = :dash,
    legend = :bottom,
    label = Int(100*gcrange[dGC[1]])
);

plot!(plot_rnafist_orffirst2, ncodons,log2.(rnafirst2./orffirst2)[dGC[2],:],
    linestyle = :solid,
    label = Int(100*gcrange[dGC[2]])
);

plot!(plot_rnafist_orffirst2, ncodons,log2.(rnafirst2./orffirst2)[dGC[3],:],
    linestyle = :dot,
    label = Int(100*gcrange[dGC[3]])
);

savefig(plot_rnafist_orffirst2, figdir*"whoisfirst2_promoterless.pdf")

# Plot single step trajectory probability versus GC content

plot_whichfirstgcr = plot(gcrange*100,log2.(rnafirst./orffirst)[:,klen],
    xlabel = "GC %",
    ylabel = "P_{ORF-first}\nP_RNA-first",
    size = (width = cm2pt(12.5), height = cm2pt(11)),
    xlim = [28,65]
);

savefig(plot_whichfirstgcr, figdir*"whoisfirstgcr_promoterless.pdf")

# Plot two step trajectory probability versus GC content

plot_whichfirstgcr2 = plot(100*gcrange,log2.(rnafirst2./orffirst2)[:,klen],
    xlabel = "GC %",
    ylabel = "P_{ORF-first}\nP_RNA-first",
    size = (width = cm2pt(12.5), height = cm2pt(11)),
    xlim = [25,62]
);

savefig(plot_whichfirstgcr2, figdir*"whoisfirstgcr2_promoterless.pdf")


plot_Loss = plot(ncodons,log2.(orfloss[dGC,:]./rnaloss[dGC,:]),
    xlabel = "ORF length (codons)",
    ylabel = "P_ORF-loss\nP_RNA-loss",
    size = (width = cm2pt(12), height = cm2pt(11)),
);

savefig(plot_Loss, figdir*"pLoss_new.pdf")

plot_loss_klen_gcrange = plot(gcrange,log2.(orfloss[:,klen]./rnaloss[:,klen]),
    xlabel = "GC-content",
    ylabel = "P_ORF-loss\nP_RNA-loss",
    size = (width = cm2pt(12.5), height = cm2pt(11)),
    ylims = [-0.3,2.1],
    xlims = [0.29,0.65]
);
savefig(plot_loss_klen_gcrange, figdir*"pLossklengcr_new.pdf")

onlyrnaloss = corfstay[dGC,:].*crnaloss[dGC,:]./corfprob[dGC,:];
onlyorfloss = crnastay[dGC,:].*corfloss[dGC,:]./crnaprob[dGC,:];

# plot_onlyLoss = plot(ncodons,log2.(onlyorfloss./onlyrnaloss),
#     xlabel = "ORF length (codons)",
#     ylabel = "P_onlyORF-loss\nP_onlyRNA-loss",
#     size = (width = cm2pt(12), height = cm2pt(11)),
# );

# savefig(plot_Loss, figdir*"pOnlyLoss_new.pdf")


# Estimate from D.mel data #

promprob = 5331240/43164203; # CRMs in open chromatin / Total intergenic region

polyavalsX = [
    sum(nprob6.(polyavars)),
    featuregain6(nopolyas,polyavars),
    featuregain6(polyavars,nopolyas)/sum(nprob6.(polyavars)),
    featurestay6(polyavars)
];

ATGvalsX = [
    nprob3("ATG"),
    featuregain3(noATG,["ATG"]),
    featuregain3(["ATG"],noATG)/nprob3("ATG"),
    featurestay3(["ATG"])
];

stopvalsX = [
    sum(nprob3.(stopvars)),
    featuregain3(nostop,stopvars),
    featuregain3(stopvars,nostop)/sum(nprob3.(stopvars)),
    featurestay3(stopvars)
];

(rnaprobX, rnagainX, rnalossX, rnastayX, 
crnaprobX, crnagainX, crnalossX, crnastayX, 
orfprobX, orfgainX, orflossX, orfstayX,
corfprobX, corfgainX, corflossX, corfstayX) = [zeros(length(ncodons),1) for i = 1:16];


for k in eachindex(ncodons)
    (rnaprobX[k], rnagainX[k], rnalossX[k], rnastayX[k]) = rnaprobs_nopromoter(polyavalsX, ncodons[k], true)

    (crnaprobX[k], crnagainX[k], crnalossX[k], crnastayX[k]) = rnaprobs_nopromoter(polyavalsX, ncodons[k], false)

    (orfprobX[k], orfgainX[k], orflossX[k], orfstayX[k]) = orfprobs(ATGvalsX,stopvalsX,ncodons[k],true,polyavalsX);

    (corfprobX[k], corfgainX[k], corflossX[k], corfstayX[k]) = orfprobs(ATGvalsX,stopvalsX,ncodons[k],false,polyavalsX);
end

rnaprobX = rnaprobX*promprob;
rnagainX = rnagainX*promprob;
crnaprobX = crnaprobX*promprob;
crnagainX = crnagainX*promprob;

genegainX = rnastayX .* corfgainX + orfstayX .* crnagainX + orfgainX.*rnagainX;

genegain2X = genegainX./(1 .-rnaprobX .- orfprobX);

genelossX = orflossX.+ rnalossX;

rnafirstX = rnastayX .* orfgainX;
orffirstX = orfstayX .* crnagainX;

onlyrnagainX = (1 .- orfprobX .- orfgainX).*rnagainX;
onlyorfgainX = (1 .- rnaprobX .- rnagainX).*orfgainX;

rnafirst2X = onlyrnagainX.*(1 .- rnalossX).*(corfgainX./(1 .-orfprobX));
orffirst2X = onlyorfgainX.*(1 .- orflossX).*(crnagainX./(1 .-rnaprobX));

plotX_genegain_geneloss = plot(ncodons,log.(genegain2X)./log.(genelossX),
    xlabel = "ORF length (codons)",
    ylabel = "# Gene Losses \n per Gene Gain",
    size = (width = cm2pt(12), height = cm2pt(11)),
    yticks = [2,2.5,3]
);
savefig(plotX_genegain_geneloss, figdir*"geneGainLoss_dmel.pdf")

plotX_rnafist_orffirst = plot(ncodons,log2.(rnafirstX./orffirstX),
    xlabel = "ORF length (codons)",
    ylabel = "P_{ORF-first}\nP_RNA-first",
    size = (width = cm2pt(12.5), height = cm2pt(11)),
    yticks = [3.5,3,2.5]
);
savefig(plotX_rnafist_orffirst, figdir*"whoisfirst_dmel.pdf")

plotX_rnafist_orffirst2 = plot(ncodons,log2.(rnafirst2X./orffirst2X),
    xlabel = "ORF length (codons)",
    ylabel = "P_{ORF-first}\nP_RNA-first",
    size = (width = cm2pt(12.5), height = cm2pt(11)),
    yticks = [0.2,0.3,0.4]
);
savefig(plotX_rnafist_orffirst2, figdir*"whoisfirst2_dmel.pdf")

plotX_Loss = plot(ncodons,log2.(orflossX./rnalossX),
    xlabel = "ORF length (codons)",
    ylabel = "P_ORF-loss\nP_RNA-loss",
    size = (width = cm2pt(12), height = cm2pt(11)),
);

savefig(plotX_Loss, figdir*"pLoss_dmel.pdf")

freegenome = readdlm("freegenome_lendist_knownchr.txt", Int32);

nsites = 4 .*[sum(freegenome[:,2].*(freegenome[:,1] .- (x + 1))) for x in (10 .+ncodons.*3)];

plot_nsites = plot(ncodons,nsites/1e7,
xlabel = "ORF length (codons)",
ylabel = "Number of feasible sites\n (x 10^7)",
size = (width = cm2pt(12), height = cm2pt(11)),
yticks = [2.6:0.4:3.8;]
);

plot_tGG = plot(ncodons,genegainX.*nsites*2,
    xlabel = "ORF length (codons)",
    ylabel = "Gene gain rate\n across multiple sites",
    size = (width = cm2pt(12.5), height = cm2pt(11)), 
    yaxis =:log, 
    yticks = [1e-12,1e-10,1e-8,1e-6]
);

savefig(plot_tGG,figdir*"pTotalGain_dmel.pdf")
savefig(plot_nsites,figdir*"pTotalSites_dmel.pdf")


# Fixation probabilities #

popsizes = 100 .* 2 .^ [1:10;];
gainfix = hcat([genegain2X/(2*x) for x in popsizes]...);
lossfix = hcat([genelossX/(2*x) for x in popsizes]...);

selcoef = 0.001

gainfixPS = hcat([genegain2X.*kimura(selcoef,2*x) for x in popsizes]...);
gainfixNS = hcat([genegain2X.*kimura(-selcoef,2*x) for x in popsizes]...);

lossfixPS = hcat([genelossX.*kimura(-selcoef,2*x) for x in popsizes]...);
lossfixNS = hcat([genelossX.*kimura(selcoef,2*x) for x in popsizes]...);

fixprobs = log.(gainfix)./log.(lossfix);
fixprobsPS = log.(gainfixPS)./log.(lossfixPS);
fixprobsNS = log.(gainfixNS)./log.(lossfixNS);


plot_fix = plot(ncodons,fixprobs[:,2:2:10],
    xlabel = "ORF length (codons)",
    ylabel = "#Gene loss with extinction \n per Gene gain with fixation",
    size = (width = cm2pt(12.5), height = cm2pt(11)),
    xlim = [25,350]
);

savefig(plot_fix, figdir*"pFix_dmel.pdf")

plot_fix = plot(ncodons,fixprobsPS[:,2:2:10],
    xlabel = "ORF length (codons)",
    ylabel = "#Gene loss with extinction \n per Gene gain with fixation",
    size = (width = cm2pt(12.5), height = cm2pt(11)),
    xlim = [25,350]
);

savefig(plot_fix, figdir*"pFix_dmel_PS.pdf")

plot_fix = plot(ncodons,fixprobsNS[:,2:2:10],
    xlabel = "ORF length (codons)",
    ylabel = "#Gene loss with extinction \n per Gene gain with fixation",
    size = (width = cm2pt(12.5), height = cm2pt(11)),
    xlim = [25,350]
);

savefig(plot_fix, figdir*"pFix_dmel_NS.pdf")

## Protein properties ##

aaprob = [sum(nucprob.(gencode[gencode[:,3] .== x,1],gccontent)) for x in aas];
aaprob2 = aaprob./sum(aaprob);
exp_hydropathy = sum(hydropathy3.*aaprob2)

acidic = indexin(["D","E"],aas);
basic = indexin(["K","R"],aas);
hydrophobic = hydropathy3 .<0;

allcodons = kmers(3);
nostopvars = allcodons[allcodons .∉ Ref(["TAA","TAG","TGA"])]

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
for i in nostopvars
    for j in nostopvars
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

