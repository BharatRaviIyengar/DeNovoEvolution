using LinearAlgebra
include("nucleotidefuncts.jl")

poisson = (l,k,n) -> l^n * exp(-l*k)/factorial(n)
poizero = (l,k) -> exp(-l*k)
geodist = (P,n) -> P*(1-P)^(n-1)

stopcodons = ["TAA","TAG","TGA"]
function nostop(codonset)
	codonset[codonset .∉ Ref(stopcodons)]
end

(xProb, xGain, xLoss, xStay) = [i for i in 1:4];

allcodons = kmers(3)
nostopcodons = nostop(allcodons);

nsub, nucsmbt = readdlm(joinpath(Base.source_dir(),organism*"_mutbias.txt"),'\t',header = true);

if organism == "scer"
    mutrate = 1.7e-10
end

if organism == "scer"
    insrate = 1.54e-12
    delrate = 3.46e-12
    ΔI = 1.875
    ΔD = 4.055

elseif organism == "dmel"
    insrate = 1.20e-10
    delrate = 4.95e-10
    ΔI = 9.083
    ΔD = 7.146
end


codonfreq = readdlm(joinpath(Base.source_dir(),organism*"_orf_codonfreq.txt"), '\t');

trimerfreq = readdlm(joinpath(Base.source_dir(),organism*"_intergenic_trimers.txt"), '\t');
trimerfreq[:,2] = normalize(trimerfreq[:,2])

noATG = setdiff(allcodons, ["ATG"]);
function ATGprobs(gccontent)
    ATGprob = nucprob("ATG",gccontent)
    ATGgain = featuregain(noATG,["ATG"],gccontent)
    ATGloss = featuregain(["ATG"],noATG,gccontent)/ATGprob
    ATGstay = featurestay(["ATG"],gccontent)
    return [ATGprob, ATGgain, ATGloss, ATGstay]
end

function ATGprobsX(t3)
    ATGprob = nprob3("ATG",t3)
    ATGgain = featuregain3(noATG,["ATG"],t3)
    ATGloss = featuregain3(["ATG"],noATG,t3)/ATGprob
    ATGstay = featurestay3(["ATG"],t3)
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

function stopprobsX(t3)
    stopprob = sum([nprob3(x,t3) for x in stopcodons])
    stopgain = featuregain3(nostopcodons,stopcodons,t3)
    stoploss = featuregain3(stopcodons,nostopcodons,t3)/stopprob
    stopstay = featurestay3(stopcodons,t3)
    return [stopprob, stopgain, stoploss, stopstay]
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

maxindelsize = 30;
fsindels = iszero.([1:maxindelsize;] .%3);
indelcodons = floor.(Int,[1:maxindelsize;]./3)
pdels = delrate.*[geodist(1/ΔD,x) for x=1:maxindelsize]
pinss = insrate.*[geodist(1/ΔI,x) for x=1:maxindelsize]
dfs = (x,y) -> fsindels[x] ? 0 : y-indelcodons[x]
ifs = (x,y) -> fsindels[x] ? indelcodons[x] : y+indelcodons[x]

function indelfun(i,j,nostp,nostopgain)
    i>j ? q = j : q = i;
    outval = 0
    for x = 1:maxindelsize
        if fsindels[x]
            if j == i+x
                outval += (i-2)*pinss[x]*(nostopgain^i)*nostp^(x/3)
            elseif i == j+x
                outval += (i-2)*pdels[x]nostopgain^(j)
            end
        else
            for k = 2:q
                outval += (pdels[x]*(nostp^(j-k-indelcodons[x]-2)) + pinss[x]*nostp^(j-k+indelcodons[x]-2))*nostopgain^(k-2)
            end
        end
    end
    return outval
end

function tprobs(i,j,ATG,stop,withindel::Bool)
    nostopstay = 1 - stop[xGain] - stop[xProb]
    nostopgain = 1 - stop[xGain]/(1-stop[xProb])
    noATGstay = 1 - ATG[xGain] - ATG[xProb]
    atgstay = 1 - ATG[xLoss]
    stopstay = 1 - stop[xLoss]
    if i==j    
        outval =  atgstay*stopstay*(nostopgain)^(i-2)
        if withindel
            outval *= poizero(insrate+delrate,i)
        end
    elseif i<j
        outval = (
            (atgstay*stop[xLoss]*stop[xStay] + stopstay*ATG[xGain])*nostopstay^(j-i-1) + 
            (j-i-1)*stop[xLoss]*stop[xStay]*ATG[xProb]*nostopstay^(j-i-2)
            )*nostopgain^(i-2)
            if withindel 
                outval += stop[xProb]*atgstay*indelfun(i,j, 1-stop[xProb],nostopgain)
            end
    else
        outval = (
            atgstay*stop[xGain]/(1-stop[xProb]) + 
            ATG[xLoss]*noATGstay^(i-j-1)*ATG[xProb]*stopstay + 
            ATG[xProb]*(i-j-1)*stop[xGain]*stopstay
            )*nostopgain^(j-2) 
            if withindel 
                outval += stop[xProb]*atgstay*indelfun(i,j, 1-stop[xProb],nostopgain)
            end
    end
    return outval
end

function T5(i,j,ATG,stop)
    nostopstay = 1 - stop[xGain] - stop[xProb]
    nostopgain = 1- stop[xGain]/(1-stop[xProb])
    noATGstay = 1 - ATG[xGain] - ATG[xProb]
    stopstay = stop[xStay]/stop[xProb]
    if i<j
        outval = (ATG[xGain]*nostopstay^(j-i-1) +(j-i-1)*stop[xLoss]*stop[xProb]*ATG[xProb]*nostopstay^(j-i-2))*nostopgain^(i-2)
    elseif i>j
        outval = (ATG[xLoss]*noATGstay^(i-j-1)*ATG[xProb] + ATG[xProb]*(i-j-1)*stop[xGain])*nostopgain^(j-2)
    else 
        outval = 0
    end
    return outval*stopstay
end

function T3(i,j,ATG,stop,withindel::Bool)
    nostopstay = 1 - stop[xGain] - stop[xProb]
    nostopgain = 1- stop[xGain]/(1-stop[xProb])
    # noATGstay = 1 - ATG[xGain] - ATG[xProb]
    atgstay = ATG[xStay]/ATG[xProb]
    if i<j
        outval = (stop[xLoss]*stop[xProb]*nostopstay^(j-i-1))*nostopgain^(i-2) 
        if withindel
            outval += stop[xProb]indelfun(i,j,stop[xProb],nostopgain)
        end
    elseif i>j
        outval = (stop[xGain]/(1-stop[xProb]))*nostopgain^(j-2)
        if withindel 
            outval += stop[xProb]*indelfun(i,j, 1-stop[xProb],nostopgain)
        end
    else 
        outval = 0
    end
    return outval*atgstay
end