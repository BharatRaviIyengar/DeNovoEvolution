using Base.Iterators
using Combinatorics
using Plots
using Measures
using LaTeXStrings

include("nucleotidefuncts.jl")

cm2pt = (cm) -> 28.3465*cm
figdir = "~/Documents/EBB/StochasticModel/Figures/"

default(linecolor = :black, linewidth = 2, tickfont = font(10), guidefontsize = 12,
        framestyle = :box, legend = false)


gccontent = 0.42

function transvals(gccontent)
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

    transprob = (tataprob + inrprob)*polyaprob
    transgain = (tatagain + inrgain)*polyastay + (tatastay + inrstay)*polyagain
    tfpinr = inrprob*(1-tataprob)/(tataprob + inrprob - inrprob*tataprob)
    tfptata = tataprob*(1-inrprob)/(tataprob + inrprob - inrprob*tataprob)
    tfpboth = inrprob*tataprob/(tataprob + inrprob - inrprob*tataprob)
    tfpboth = 1 - tfptata - tfpinr

    transloss = ( tfpinr*inrloss2*notatastay +
      tfptata*tataloss2*noinrstay +
      tfpboth*tataloss2*inrloss2 +
      polyaloss2
    )
    transstay = (inrstay + tatastay)*polyastay
    return [transprob; transgain; transloss; transstay]
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


ncodons = [40:300;]

ATGvalues = ATGvals(gccontent)
stopvalues = stopvals(gccontent)

(orfprob, orfgain, orfloss, orfstay) = [zeros(size(ncodons)) for i = 1:4 ]

for k in 1:length(ncodons)
    (orfprob[k], orfgain[k], orfloss[k], orfstay[k]) = orfvals(gccontent,ATGvalues,stopvalues,ncodons[k])
end

genegain = (transgain + transstay) .* orfgain + orfstay .* transgain
geneloss = orfloss + transloss

rnafirst = transstay .* orfgain
orffirst = orfstay .* transgain

plot_genegain_geneloss = plot(ncodons,log.(genegain)./log.(geneloss),
    xlabel = "ORF length (codons)",
    ylabel = "Gene losses per gene gain",
    size = (width = cm2pt(30), height = cm2pt(11)),
)
savefig(pComb, figdir*"TvsO.pdf")

tl = log10.(transloss./orfgain)
ol = log10.(orfloss/transgain)

# pP = plot(ncodons,transprob./orfprob,
#     yaxis=:log,
#     title = "Stationary probability",
#     xlabel = "ORF length",
#     ylabel = L"P_\textrm{RNA}/P_\textrm{ORF}",
#     xrotation = 45,
# );

# pG = plot(ncodons,transgain./orfgain,
#     yaxis=:log,
#     title = "Gain probability",
#     xlabel = "ORF length",
#     ylabel = L"P_\textrm{RNA}/P_\textrm{ORF}",
#     xrotation = 45
# );

# pL = plot(ncodons,transloss./orfloss,
#     title = "Loss probability",
#     xlabel = "ORF length",
#     ylabel = L"P_\textrm{RNA}/P_\textrm{ORF}",
#     xrotation = 45
# );
# pComb = plot!(pP,pG,pL,
#     size = (width = cm2pt(30), height = cm2pt(11)),
#     left_margin = 3.5mm,
#     bottom_margin = 7mm,
#     layout = @layout [a b c]
#     );

# savefig(pComb, figdir*"TvsO.pdf")
