using Base.Iterators
using Combinatorics
using Plots
using Measures
using LaTeXStrings

include("nucleotidefuncts.jl")

cm2pt = (cm) -> 28.3465*cm
figdir = joinpath(Base.source_dir(),"../Manuscripts/Figures/");

default(linecolor = :black, linewidth = 2, tickfont = font(10,"Helvetica"), 
guidefont = font(13,"Helvetica"),framestyle = :box, legend = false);


gccontent = 0.42

function rnavals(gccontent)
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
    ))

    notata = setdiff(kmers(8),tatavars)

    tataprob = sum([nucprob(x,gccontent) for x in tatavars])
    tatagain = featuregain(notata,tatavars,gccontent)
    tataloss1 = featuregain(tatavars,notata,gccontent)
    tataloss2 = tataloss1/tataprob
    tatastay = featurestay(tatavars,gccontent)
    notatastay = 1 - tatagain - tataprob

    ## Initiator element
    # Consensus sequence = "BBCABW";  DOI:10.1101/gad.295980.117

    mers6 = kmers(6)
    inrvars = vec(join.(product("TGC","TGC","C", "A","TGC","TA")))
    noinrs = setdiff(mers6, inrvars)
    inrprob = sum([nucprob(x,gccontent) for x in inrvars])
    inrgain = featuregain(noinrs,inrvars,gccontent)
    inrloss1 = featuregain(inrvars,noinrs,gccontent)
    inrloss2 = inrloss1/inrprob
    inrstay = featurestay(inrvars,gccontent)
    noinrstay = 1- inrgain -inrprob

    ## PolyA signal

    polyavars = ["AATAAA";"ATTAAA"; "AGTAAA";"TATAAA"]
    nopolyas = setdiff(mers6, polyavars)
    polyaprob = sum([nucprob(x,gccontent) for x in polyavars])
    polyagain = featuregain(nopolyas,polyavars,gccontent)
    polyaloss1 = featuregain(polyavars,nopolyas,gccontent)
    polyaloss2 = polyaloss1/polyaprob
    polyastay = featurestay(polyavars,gccontent)

    ## Transcription
    # (Initiator or TATA) and polyA

    rnaprob = (tataprob + inrprob)*polyaprob
    rnagain = (tatagain + inrgain)*polyastay + (tatastay + inrstay)*polyagain
    tfpinr = inrprob*(1-tataprob)/(tataprob + inrprob - inrprob*tataprob)
    tfptata = tataprob*(1-inrprob)/(tataprob + inrprob - inrprob*tataprob)
    tfpboth = inrprob*tataprob/(tataprob + inrprob - inrprob*tataprob)
    tfpboth = 1 - tfptata - tfpinr

    rnaloss = ( tfpinr*inrloss2*notatastay +
      tfptata*tataloss2*noinrstay +
      tfpboth*tataloss2*inrloss2 +
      polyaloss2
    )
    rnastay = (inrstay + tatastay)*polyastay
    return [rnaprob; rnagain; rnaloss; rnastay]
end

function ATGvals(gccontent)
    ## Start codon
    noATG = setdiff(kmers(3), ["ATG"])
    ATGprob = nucprob("ATG",gccontent)
    ATGgain = featuregain(noATG,["ATG"],gccontent)
    ATGloss1 = featuregain(["ATG"],noATG,gccontent)
    ATGloss2 = ATGloss1/ATGprob
    ATGstay = featurestay(["ATG"],gccontent)
    return [ATGprob, ATGgain, ATGloss1, ATGloss2, ATGstay]
end

function stopvals(gccontent)
    ## Stop codon
    stopvars = ["TAA","TAG","TGA"]
    nostop = setdiff(kmers(3), stopvars)
    stopprob = sum([nucprob(x,gccontent) for x in stopvars])
    stopgain = featuregain(nostop,stopvars,gccontent)
    stoploss1 = featuregain(stopvars,nostop,gccontent)
    stoploss2 = stoploss1/stopprob
    stopstay = featurestay(stopvars,gccontent)
    nostopstay = featurestay(nostop,gccontent)
    return [stopprob, stopgain, stoploss1, stoploss2, stopstay, nostopstay]
end



function orfvals(gccontent,ATGvalues,stopvalues,k)
    ## Open reading frame
    (ATGprob,ATGgain, ATGloss2, ATGstay) = ATGvalues[[1,2,4,5]]
    (stopprob, stopgain, stoploss1, stoploss2, stopstay, nostopstay) = stopvalues
    orfprob = ATGprob*stopprob*(1 - stopprob)^(k-2)
    orfgain = nostopstay^(k-2)*(
stopstay*ATGgain + ATGstay*stopgain +
(k-2)*(stopstay*ATGstay*stoploss1/nostopstay)
)
    orfloss = stoploss2 + ATGloss2 + (k-2)*stopgain/(1-stopprob)
    orfstay = ATGstay*stopstay*(nostopstay)^(k-2)
    return [orfprob; orfgain; orfloss; orfstay]
end


## Protogene
# Has a length of x
# Includes promoters and polyA signal
# Can produce a RNA of length y ≤ x – 2 – 6 (–32; TATA)
# Can have ORF of length 3k ≤ x (IGORF) or 3k ≤ y (protein coding)

# z = 10000
# nmer = (z,x) -> z - x + 1
# yi = (x) -> x - 8
# yt = (x) -> x - 40
# nORF = (y,k) -> y - 4 - 3*k + 1
# minRNAlen = (k) -> 3*k + 4


ncodons = [30:300;];

ATGvalues = ATGvals(gccontent);
stopvalues = stopvals(gccontent);


(rnaprob, rnagain, rnaloss, rnastay) = rnavals(gccontent)
(orfprob, orfgain, orfloss, orfstay) = [zeros(size(ncodons)) for i = 1:4 ]

for k in eachindex(ncodons)
    (orfprob[k], orfgain[k], orfloss[k], orfstay[k]) = orfvals(gccontent,ATGvalues,stopvalues,ncodons[k])
end

genegain = ((rnagain + rnastay) .* orfgain + orfstay .* rnagain);
genegain2 = genegain./(1-rnaprob .- orfprob);
geneloss = orfloss .+ rnaloss;

rnafirst = rnastay .* orfgain;
orffirst = orfstay .* rnagain;

onlyrnagain = (1 .- orfprob .- orfgain).*rnagain;
onlyorfgain = (1 -rnaprob - rnagain).*orfgain;

rnafirst2 = onlyrnagain*(1 - rnaloss).*(orfgain./(1 .-orfprob));
orffirst2 = onlyorfgain.*(1 .-orfloss).*(rnagain/(1-rnaprob));

plot_genegain_geneloss = plot(ncodons,log.(genegain2)./log.(geneloss),
    xlabel = "ORF length (codons)",
    ylabel = "# Gene Losses \n per Gene Gain",
    size = (width = cm2pt(12), height = cm2pt(11)),
);
savefig(plot_genegain_geneloss, figdir*"geneGainLoss.pdf")

plot_rnafist_orffirst = plot(ncodons[1:50],log.(orffirst./rnafirst)[1:50],
    xlabel = "ORF length (codons)",
    ylabel = "P_{ORF-first}\nP_RNA-first",
    size = (width = cm2pt(12), height = cm2pt(11)),
    xticks = [35:7:80;]
);
savefig(plot_rnafist_orffirst, figdir*"first_ORF_RNA.pdf")

plot_rnafist_orffirst2 = plot(ncodons[1:50],log.(orffirst2./rnafirst2)[1:50],
    xlabel = "ORF length (codons)",
    ylabel = "P_{ORF-first}\nP_RNA-first",
    size = (width = cm2pt(12), height = cm2pt(11)),
    xticks = [35:7:80;]
);
savefig(plot_rnafist_orffirst2, figdir*"first_ORF_RNA2.pdf")

tl = log10.(rnaloss./orfgain);
ol = log10.(orfloss/rnagain);

plot_Loss = plot(ncodons[1:85],log.(orfloss./rnaloss)[1:85],
    xlabel = "ORF length (codons)",
    ylabel = "P_RNA-loss\nP_ORF-loss",
    size = (width = cm2pt(12), height = cm2pt(11)),
    xticks = [33:11:110;]
);

savefig(plot_Loss, figdir*"pLoss.pdf")

onlyrnaloss = orfstay .*(rnaprob*rnaloss)

aaprob = [sum(nucprob.(gencode[gencode[:,3] .== x,1],gccontent)) for x in aas];
aaprob2 = aaprob./sum(aaprob);
exp_hydropathy = sum(hydropathy.*aaprob2)

acidic = indexin(["D","E"],aas);
basic = indexin(["K","R"],aas);
hydrophobic = hydropathy .>0;

pHydrophobic = sum(aaprob[hydrphobic])

plot_aafreq = bar(aas,aaprob2,
    xticks = (0.5:20, aas),
    fill = :black,
    xlabel = "Amino acid",
    ylabel = "Expected frequency",
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

aasubprob_conditional = aasubprob./aaprob2

asymmetric_transitions = log2.(aasubprob./aasubprob');

asymmetric_conditional = log2.(aasubprob_conditional./aasubprob_conditional');

toHydrophobic = sum(aasubprob[:,hydrophobic]);
frHydrophobic = sum(aasubprob[hydrophobic,:]);

toHydrophobic_cond = sum(aasubprob_conditional[:,hydrophobic]);
frHydrophobic_cond = sum(aasubprob_conditional[hydrophobic,:]);

plot_aasubprob = heatmap(aas,aas,log10.(aasubprob./aaprob), 
    xticks = (0.5:20, aas), 
    yticks = (0.5:20,aas),
    size = (width = cm2pt(15), height = cm2pt(14)),
    colorbar_title = "Substitution Probability",
    legend = :bottom
);

savefig(plot_aasubprob, figdir*"pAASub.pdf")

plot_aasymprob = heatmap(aas,aas,asymmetric_transitions, 
    xticks = (0.5:20, aas), 
    yticks = (0.5:20,aas),
    size = (width = cm2pt(15), height = cm2pt(14)),
    colorbar_title = "Substitution likelihood",
    legend = :bottom,
    yflip = true,
    ylabel = "Original amino acid",
    xlabel = "Substituted amino acid",
    c = :bluesreds
);

savefig(plot_aasymprob, figdir*"pAAsubSym.pdf")

plot_aasymprob2 = heatmap(aas,aas,asymmetric_conditional, 
    xticks = (0.5:20, aas), 
    yticks = (0.5:20,aas),
    size = (width = cm2pt(15), height = cm2pt(14)),
    colorbar_title = "Substitution likelihood",
    legend = :bottom,
    yflip = true,
    ylabel = "Original amino acid",
    xlabel = "Substituted amino acid",
    c = :bluesreds
);


savefig(plot_aasymprob2, figdir*"pAAsubSymCond.pdf")