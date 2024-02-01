using StatsBase
include("nucleotidefuncts.jl")

organism = "scer"

nsub, nucsmbt = readdlm(joinpath(Base.source_dir(),organism*"_mutbias.txt"),'\t',header = true);

if(organism=="scer")
    promprob = 0.458234
    mutrate = 1.7e-10
elseif(organism =="dmel")
    promprob = 5331240/43164203; # CRMs in open chromatin / Total intergenic region
    mutrate = 7.8e-9
end

freegenome = readdlm(joinpath(Base.source_dir(),organism*"_intergenic_lendist.txt"),'\t',Int);

nsites = (x) -> sum(max.(0,freegenome[:,1] .-x .+ 1).*freegenome[:,2])

allcodons = kmers(3);

noATG = setdiff(allcodons, ["ATG"]);
stopvars = ["TAA","TAG","TGA"];
nostop = setdiff(allcodons, stopvars);
polyavars = ["AATAAA";"ATTAAA"; "AGTAAA";"TATAAA"];
nopolyas = setdiff(kmers(6), polyavars);

f6 = readdlm(joinpath(Base.source_dir(),organism*"_intergenic_hexamers.txt"), '\t');
f3 = readdlm(joinpath(Base.source_dir(),organism*"_intergenic_trimers.txt"), '\t');

(xProb, xGain, xLoss, xStay) = [i for i in 1:4];

function rnaprobs(promprob,polya,rnalen)
    
    noPAsites = rnalen - 5 -6 # total 6-mers in the interior  

    nopolyaprob = 1 - polya[xProb]
    nopolyastay = nopolyaprob  - polya[xGain]

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

    rnagain =  promprob*(gainmech1 + gainmech2) 
    
    rnaloss = polya[xLoss] + noPAsites*polya[xGain];

    rnastay = promprob*polya[xStay]*nopolyastay^noPAsites;

    return [rnaprob; rnagain; rnaloss; rnastay]
end

function tprobs(i,j,ATG,stop)
    nostopstay = 1 - stop[xGain] - stop[xProb]
    noATGstay = 1 - ATG[xGain] - ATG[xProb]
    if i==j    
        return (1-ATG[xLoss])*(1-stop[xLoss])*(1-stop[xGain]/stop[xProb])^i
    elseif i<j
        return (stop[xLoss] + ATG[xGain])*nostopstay^(j-i)
    else
        return stop[xGain]/stop[xProb] + ATG[xLoss]*noATGstay^(i-j)*ATG[xProb]
    end
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



function polyaprobs(f6)

    ppa = sum([nprob6(x,f6) for x in polyavars]);

    # u1npa = vcat(overlapXnucs.(nopolyas,-1)...); # nucleotides 1nt upstream of non-polya

    upa = [vcat(overlapXnucs.(polyavars,-x)...) for x = 1:5]; # nucleotides upstream of a polya

    ppa_upa = zeros(5)
    gpa_upa = zeros(5)
    for y=1:5
        ixx = upa[y] .∈ Ref(polyavars)
        if(!any(ixx))
            ppa_upa[y] = 0
        else
        ppa_upa[y] = sum(nprob6(x,f6) for x in upa[y][ixx])/sum(nprob6(x,f6) for x in upa[y]);
        end
        gpa_upa[y] = featuregain6(upa[y][.!ixx],polyavars,f6);
    end
    gpa = featuregain6(nopolyas,polyavars,f6);
    lpa = featuregain6(polyavars,nopolyas,f6)/sum([nprob6(x,f6) for x in polyavars]);
    spa = featurestay6(polyavars,f6);

    basic = [ppa,gpa,lpa,spa]

    return basic, ppa_upa, gpa_upa

end

function rnaprob(promprob,rnalen,ppx_bundle)
    prob_basic = ppx_bundle[1]
    ppa_overlaps = ppx_bundle[2]
    internal_sites = rnalen-12+1
    
    rnaprob = polya[xProb] * nopolyaprob^noPAsites;

    prob_basic[xProb] * (1-prob_basic[xProb])^internal_sites * prod(1-ppa_overlaps[x] for x=1:5)
end



ATGvalsX = [
    nprob3("ATG",f3),
    featuregain3(noATG,["ATG"],f3),
    featuregain3(["ATG"],noATG,f3)/nprob3("ATG",f3),
    featurestay3(["ATG"],f3)
];

stopvalsX = [
    sum([nprob3(x,f3) for x in stopvars]),
    featuregain3(nostop,stopvars,f3),
    featuregain3(stopvars,nostop,f3)/sum([nprob3(x,f3) for x in stopvars]),
    featurestay3(stopvars,f3)
];

ncodons = [30:300;]; # max > 99.9999 % CI
lrangeRNA = [100:2000;]; # max > 99.9999 % CI

(rnaprobX, rnagainX, rnalossX, rnastayX) = [zeros(length(lrangeRNA),1) for i = 1:4]

(orfprobX, orfgainX, orflossX, orfstayX,
corfprobX, corfgainX, corflossX, corfstayX) = [zeros(length(ncodons),1) for i = 1:8];


for k in eachindex(lrangeRNA)
    (rnaprobX[k], rnagainX[k], rnalossX[k], rnastayX[k]) = 
    rnaprobs(promprob,polyavalsX,k);
end

for k in eachindex(ncodons)
    (orfprobX[k], orfgainX[k], orflossX[k], orfstayX[k]) = orfprobs(ATGvalsX,stopvalsX,ncodons[k],true,polyavalsX);

    (corfprobX[k], corfgainX[k], corflossX[k], corfstayX[k]) = orfprobs(ATGvalsX,stopvalsX,ncodons[k],false,polyavalsX);
end

## Expected rate that at least one RNA is gained anywhere in the intergenic region (of different lengths), or at least one ORF is gained in any intergenic RNA ## 

no_rna_gains = 1 .- rnagainX;

p_no_orfs= (l) -> prod([(1-corfprobX[x])^(l-ncodons[x]+1) for x in findall(ncodons .< l/3)])

g_no_orfs = (l) -> prod([(1-corfgainX[x])^(l-ncodons[x]+1) for x in findall(ncodons .< l/3)])

any_orf_gain_in_RNA = 1 .- rnaprobX.*[g_no_orfs(x) for x in lrangeRNA];

# rate of any such event (differnt lengths)
mu = (rnaprobX.*(1 .-no_orf_gain_in_RNA .*rnalossX) .+ rnagainX).^nsites.(lrangeRNA);

## Probability at least one such event happens ##
px = 1- prod(exp.(-mu)) 

px2 = px +  0.5*(1-exp(-mutrate*sum(ncodons.*rnaprobX.*corfprobX.*nsites[:,2]*2)))


function runsim(popsize,generations)
    rn_muts = rand(popsize,generations)

end

