using LinearAlgebra
using Measures
using Plots
using Printf

cd(Base.source_dir())
include("nucleotidefuncts.jl")

organism = "dmel"
# organism = "scer"

poisson = (l,k,n) -> l^n * exp(-l*k)/factorial(n)


cm2pt = (cm) -> 28.3465*cm
figdir = joinpath(Base.source_dir(),"../Manuscripts/M3/Figures/");
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

pnomut = (l,m) -> exp(-l*m)

normalize = x -> x/sum(x)

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

function tprobs(i,j,ATG,stop)
    nostopstay = 1 - stop[xGain] - stop[xProb]
    noATGstay = 1 - ATG[xGain] - ATG[xProb]
    if i==j    
        return (1-ATG[xLoss])*(1-stop[xLoss])*(1-stop[xGain]/stop[xProb])^i
    elseif i<j
        return (stop[xLoss]*stop[xProb] + ATG[xGain])*nostopstay^(j-i) +(j-i)*stop[xLoss]*stop[xProb]*ATG[xProb]*nostopstay^(j-i-1)
    else
        return stop[xGain]/(1-stop[xProb]) + ATG[xLoss]*noATGstay^(i-j)*ATG[xProb] + ATG[xProb]*(i-j)*stop[xGain]
    end
end

function T5(i,j,ATG,stop)
    nostopstay = 1 - stop[xGain] - stop[xProb]
    noATGstay = 1 - ATG[xGain] - ATG[xProb]
    if i<j
        return ATG[xGain]*nostopstay^(j-i) +(j-i)*stop[xLoss]*stop[xProb]*ATG[xProb]*nostopstay^(j-i-1)
    elseif i>j
        return ATG[xLoss]*noATGstay^(i-j)*ATG[xProb] + ATG[xProb]*(i-j)*stop[xGain]
    else return 1
    end
end

function T3(i,j,ATG,stop)
    nostopstay = 1 - stop[xGain] - stop[xProb]
    noATGstay = 1 - ATG[xGain] - ATG[xProb]
    if i<j
        return stop[xLoss]*stop[xProb]*nostopstay^(j-i)
    elseif i>j
        return stop[xGain]/(1-stop[xProb])
    else return 1
    end
end

function markovcalc(gens,transmat,orfvalues,traj="s")
    xmt = transmat^gens;
    if traj=="e"
        xinit = orfvalues[xProb,:]'.*LowerTriangular(ones(size(transmat)));
        setindex!.(Ref(xinit), 0.0, eachindex(ncodons), eachindex(ncodons));
    elseif traj=="t"
        xinit = UpperTriangular(ones(size(transmat))).*orfvalues[xProb,:]';
        setindex!.(Ref(xinit), 0.0, eachindex(ncodons), eachindex(ncodons));
    else
        xinit = orfvalues[xProb,:].*Matrix(1.0*I,size(transmat));
    end
    return xinit*xmt
end

gcrange = [0.3:0.1:0.6;];
#ncodons = [20:500;];

ncodons =[3:900;];

similarity90 = [pnomut(mutrate,x) for x in Int.(ceil.(ncodons.*0.9))];

atgvals = hcat([ATGprobs(g) for g in gcrange]...);
stopvals = hcat([stopprobs(g) for g in gcrange]...);

orfvals = zeros(4,length(ncodons),length(gcrange));
tmat = zeros(length(ncodons),length(ncodons),length(gcrange))

m3 = zeros(length(ncodons),length(ncodons),length(gcrange))
m5 = zeros(length(ncodons),length(ncodons),length(gcrange))

for g in eachindex(gcrange)
    orfvals[:,:,g] = hcat([orfprobs(atgvals[:,g],stopvals[:,g],k) for k in ncodons]...);
    tmat[:,:,g] = similarity90.*[tprobs(i,j,atgvals[:,g],stopvals[:,g]) for i in ncodons, j in ncodons];
    m5[:,:,g] = [T5(i,j,atgvals[:,g],stopvals[:,g]) for i in ncodons, j in ncodons]
    m3[:,:,g] = [T3(i,j,atgvals[:,g],stopvals[:,g]) for i in ncodons, j in ncodons]
end
atgvalsX = ATGprobsX(trimerfreq);
stopvalsX = stopprobsX(trimerfreq);
orfvalsX = hcat([orfprobs(atgvalsX,stopvalsX,k) for k in ncodons]...);

tmatX = similarity90.*[tprobs(i,j,atgvalsX,stopvalsX) for i in ncodons, j in ncodons];

m3x = [T3(i,j,atgvalsX,stopvalsX) for i in ncodons, j in ncodons];
m5x = [T5(i,j,atgvalsX,stopvalsX) for i in ncodons, j in ncodons];

m5g3x = log.(m5x./m3x);

# m5g3x[findall(m5g3x.<=0)] .= NaN;

mostlikelystart = zeros(Int32,length(ncodons),length(gcrange)+1);
for g in eachindex(gcrange)
    egt = eigen(tmat[:,:,g]);
    vecs = orfvals[xProb,:,g].*egt.vectors; 
    vecs2 = abs.(vecs)
    mostlikelystart[:,g] = [ncodons[findlast(vecs2[:,x] .== maximum(vecs2[:,x]))] for x in eachindex(ncodons)];
end

egt = eigen(tmatX);
vecs = orfvalsX[xProb,:].*egt.vectors; 
vecs2 = abs.(vecs)
mostlikelystart[:,end] = [ncodons[findlast(vecs2[:,x] .== maximum(vecs2[:,x]))] for x in eachindex(ncodons)];

gens = 2 .^[0:Int(floor(log2(10/mutrate)));]
(EgT, EgS, TgS) = [zeros(Int64,length(gens),length(gcrange)+1) for i=1:3];

(ev, tv, sv, Eop, Top, Sop) = [zeros(length(ncodons),length(gens),length(gcrange)+1) for i=1:6];

for r in eachindex(gens)
    for g in eachindex(gcrange)
        tmg = tmat[:,:,g]^gens[r]
        sameQ = markovcalc(1,tmg, orfvals[:,:,g],"s");
        extnQ = markovcalc(1,tmg, orfvals[:,:,g],"e");
        trncQ = markovcalc(1,tmg, orfvals[:,:,g],"t");
        (sv[:,r,g], ev[:,r,g], tv[:,r,g]) = [diag(x) for x in [sameQ, extnQ, trncQ]];
        # println([gcrange[g],gens[r]])
        if(any(ev[:,r,g] .> tv[:,r,g]))
            EgT[r,g] = ncodons[findfirst(ev[:,r,g] .> tv[:,r,g])]
        end    
        if(any(ev[:,r,g] .> sv[:,r,g]))
            EgS[r,g] = ncodons[findfirst(ev[:,r,g] .> sv[:,r,g])]
        end 
        if(any(tv[:,r,g] .> sv[:,r,g]))
            TgS[r,g] = ncodons[findfirst(tv[:,r,g] .> sv[:,r,g])]
        end 
        Eop[:,r,g] = [sum(tmg[x,x+1:end]) for x in eachindex(ncodons)]
        Top[:,r,g] = [sum(tmg[x,1:x-1]) for x in eachindex(ncodons)]
        Sop[:,r,g] = diag(tmg)
    end
    tmg = tmatX^gens[r]
    sameQ = markovcalc(1,tmg,orfvalsX,"s");
    extnQ = markovcalc(1,tmg,orfvalsX,"e");
    trncQ = markovcalc(1,tmg,orfvalsX,"t");
    (sv[:,r,end], ev[:,r,end], tv[:,r,end]) = [diag(x) for x in [sameQ, extnQ, trncQ]];
    if(any(ev[:,r,end] .> tv[:,r,end]))
        EgT[r,end] = ncodons[findfirst(ev[:,r,end] .> tv[:,r,end])]
    end    
    if(any(ev[:,r,end] .> sv[:,r,end]))
        EgS[r,end] = ncodons[findfirst(ev[:,r,end] .> sv[:,r,end])]
    end
    if(any(tv[:,r,end] .> sv[:,r,end]))
        TgS[r,end] = ncodons[findfirst(tv[:,r,end] .> sv[:,r,end])]
    end
    Eop[:,r,end] = [sum(tmg[x,x+1:end]) for x in eachindex(ncodons)]
    Top[:,r,end] = [sum(tmg[x,1:x-1]) for x in eachindex(ncodons)]
    Sop[:,r,end] = diag(tmg)
    println(r)
end


fac = Int(5*10^(Int(floor(-log10(mutrate))) -2))
eggen = findfirst(gens .>= fac)

# sameQ = markovcalc(fac,tmatX, orfvalsX,"s");
# extnQ = markovcalc(fac,tmatX, orfvalsX,"e");
# trncQ = markovcalc(fac,tmatX, orfvalsX,"t");
# (sv, ev, tv) = [diag(x) for x in [sameQ, extnQ, trncQ]];


# Plots #

pEV = plot( 
    xlabel = "Final ORF length (codons)",
    ylabel = "Initial ORF length (codons)",
    size = (width = cm2pt(11), height = cm2pt(10))
    );

for q in 1:5
    plot!(pEV,ncodons,mostlikelystart[:,q],
        linecolor = colors[q]
    );
end

mx = round(maximum(m5g3x),digits=2);
mn = round(minimum(m5g3x),digits=2);
rrx = mx - mn;
cr = cgrad(["#d40000",:black,"#5599ff"], [-mn/rrx, 1])
setindex!.(Ref(m5g3x), NaN, eachindex(ncodons), eachindex(ncodons));

plt5g3x = heatmap(m5g3x, 
    c = cr,
    clim = (mn, mx),
    cb = :top,
    xrotation = 90,
    xticks = [100:100:900;],
    yticks = [100:100:900;],
    ylabel = "Initial ORF length",
    xlabel = "Final ORF length",
    size = (width = cm2pt(11.5), height = cm2pt(11)),
);

savefig(plt5g3x,figdir*"5g3_"*organism*".pdf")


TgE_gORF = hcat([ncodons[[findfirst(Top[:,x,g] .> Eop[:,x,g]) for x in eachindex(gens)]] for g =1:5]...);

pltTgE = plot( 
    xlabel = "Generations",
    ylabel = "Minimum ORF length\n(codons)",
    title = "P(T) > P(E)",
    xaxis = :log10,
    xticks = 10 .^[0:2:10;],
    size = (width = cm2pt(11.5), height = cm2pt(11))
    );

for q in 1:5
    plot!(pltTgE,gens,TgE_gORF[:,q],
        linecolor = colors[q]
    );
end

savefig(pltTgE,figdir*"TgEvGen_fromORF_"*organism*".pdf")

pltEgT = plot( 
    xlabel = "Generations",
    ylabel = "Minimum ORF length\n(codons)",
    title = "P(E) > P(T)",
    xaxis = :log10,
    xticks = 10 .^[0:2:10;],
    size = (width = cm2pt(11.5), height = cm2pt(11))
    );

for q in 1:5
    plot!(pltEgT,gens,EgT[:,q],
        linecolor = colors[q]
    );
end

savefig(pltEgT,figdir*"EgTvGen_"*organism*".pdf")

# SgE_gORF = hcat([ncodons[[findfirst(Sop[:,x,g] .> Eop[:,x,g]) for x in eachindex(gens)]] for g =1:5]...);

# pltSgE = plot( 
#     xlabel = "Generations",
#     ylabel = "Minimum ORF length\n(codons)",
#     title = "P(S) > P(E)",
#     xaxis = :log10,
#     xticks = 10 .^[0:2:10;],
#     size = (width = cm2pt(11.5), height = cm2pt(11))
#     );

# for q in 1:5
#     plot!(pltSgE,gens,SgE_gORF[:,q],
#         linecolor = colors[q]
#     );
# end

# savefig(pltTgE,figdir*"SgEvGen_fromORF_"*organism*".pdf")


mESGx = [findfirst(EgS[:,x] .>0) for x = 1:5];
mTSGx = [findfirst(TgS[:,x] .>0) for x = 1:5];

pltEgS = plot( 
    xlabel = "Generations",
    ylabel = "Minimum ORF length\n(codons)",
    title = "P(E) > P(C)",
    xaxis = :log10,
    xticks = 10 .^[floor(log10(gens[minimum(mESGx)])):10;],
    size = (width = cm2pt(11.5), height = cm2pt(11))
    );

for q in 1:5
    plot!(pltEgS,gens[mESGx[q]:end],EgS[mESGx[q]:end,q],
        linecolor = colors[q]
    );
end

savefig(pltEgS,figdir*"EgSvGen_"*organism*".pdf")


pltTgS = plot( 
    xlabel = "Generations",
    ylabel = "Minimum ORF length\n(codons)",
    title = "P(T) > P(C)",
    xaxis = :log10,
    xticks = 10 .^[floor(log10(gens[minimum(mTSGx)])):10;],
    size = (width = cm2pt(11.5), height = cm2pt(11))
    );

for q in 1:5
    plot!(pltTgS,gens[mTSGx[q]:end],TgS[mTSGx[q]:end,q],
        linecolor = colors[q]
    );
end

savefig(pltTgS,figdir*"TgSvGen_"*organism*".pdf")

plt_eg = plot(ncodons[1:300],log10.(sv)[1:300], linecolor = :black, label = "Constant", xlabel = "ORF length (codons)", ylabel = "log10(Probability)")
plot!(plt_eg, ncodons[1:300],log10.(ev)[1:300], linecolor = :red, label = "Extended");
plot!(plt_eg,ncodons[1:300],log10.(tv)[1:300], linecolor = :blue, label = "Truncated");

plot!(plt_eg, legend = :bottom, legendfont = font(11,"Helvetica"),
    size = (width = cm2pt(13), height = cm2pt(11))
    )

savefig(plt_eg,figdir*organism*"_Prob_"* (@sprintf "%0.0e" fac)*"gens.pdf")

allplots = []

for i in Int.(floor.([1:length(gens)/4:length(gens);]))
    plt_egOP = plot(ncodons[2:300],log10.(Eop[2:300,i,5]), linecolor = :red, title = string(gens[i])*" generations", label ="Extended");

    plot!(plt_egOP,ncodons[2:300],log10.(Top[2:300,i,5]), linecolor = :blue, label ="Truncated")

    if i==27
        plot!(plt_egOP, legend = :bottom, legendfont = font(11,"Helvetica"),
        size = (width = cm2pt(13), height = cm2pt(11))
        )
    end
    push!(allplots,plt_egOP)
end

apx = plot(allplots..., layout = (2,2),size = (width = cm2pt(25), height = cm2pt(20)));

savefig(apx,figdir*"Supp/TvE_fromORF"*organism*".pdf")


plots_5vs3 = Array{Plots.Plot{Plots.GRBackend}}(undef,1,4);

m35 = log2.(m3./m5);

m35[findall(m35.<=0)] .= NaN;

for i in eachindex(gcrange)
    plots_5vs3[i] = heatmap(m35[:,:,i], legend = :bottom, 
    xrotation = 90,
    ylabel = "Initial ORF length",
    xlabel = "Final ORF length",
    title = "GC% = "*string(Int(gcrange[i]*100))
    );
end

P53 = plot(plots_5vs3..., layout = (2,2), size = (width = cm2pt(25), height = cm2pt(20)));

savefig(P53,figdir*"Supp/ET_3v5_"*organism*".pdf")
