using Plots
using Measures
using LinearAlgebra
using Printf
# using LaTeXStrings
cd(Base.source_dir())
include("nucleotidefuncts.jl")

organism = "dmel"

cm2pt = (cm) -> 28.3465*cm
figdir = joinpath(Base.source_dir(),"../Manuscripts/Other_Figures/GeneLength/");
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

using DelimitedFiles

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
    ORFv = orfprobs(ATG,stop,i)
    if i==j    
        return (1-ATG[xLoss])*(1-stop[xLoss])*(1-stop[xGain]/stop[xProb])^i
    elseif i<j
        return (stop[xLoss] + ATG[xGain])*nostopstay^(j-i)
    else
        return stop[xGain]/stop[xProb] + ATG[xLoss]*noATGstay^(i-j)*ATG[xProb]
    end
end


 ncodons = [20:500;];
# gcrange = [0.3:0.05:0.6;]
# crossover = zeros(Int64,length(ncodons),length(gcrange),2);
# crossoverX = zeros(Int64,length(ncodons),2);
atgvalsX = ATGprobsX(trimerfreq);
stopvalsX = stopprobsX(trimerfreq);
orfvalsX = hcat([orfprobs(atgvalsX,stopvalsX,k) for k in ncodons]...);
tmatX = [tprobs(i,j,atgvalsX,stopvalsX) for i in ncodons, j in ncodons];
# tmatX2 = tmatX/sum(tmatX);

egt = eigen(tmatX);
vecs = orfvalsX[xProb].*egt.vectors; 
vecs2 = abs.(vecs)
mostlikelystart = [ncodons[findlast(vecs2[:,x] .== maximum(vecs2[:,x]))] for x in eachindex(ncodons)];

px = plot(ncodons[1:end-1],mostlikelystart[1:end-1], 
    xlabel = "Final ORF length (codons)",
    ylabel = "Initial ORF length (codons)"
    );

eq = ncodons[findfirst(mostlikelystart .<= ncodons)];
scatter!(px, [eq],[eq])

# for x in eachindex(ncodons)

#     (xinit1,xinit2) = [orfvalsX[xProb,:] for i=1:2];
#     xinit1[x:end] .= 0;
#     xinit2[1:(x-1)] .= 0;
#     x1 = [tmatX[x,x]^n *orfvalsX[xProb,x] for n=80:120];
#     x2 = vcat([(xinit1' * tmatX^n)[x] for n=80:120]...);
#     x3 = vcat([(xinit2' * tmatX^n)[x] for n=80:120]...);

#     c1 = findfirst(x2./x1 .>= 1.1)
#     c2 = findfirst(x3./x1 .>= 1.1)

#     isnothing(c1) ? crossoverX[x,1] = 0 : crossoverX[x,1] = c1
#     isnothing(c2) ? crossoverX[x,2] = 0 : crossoverX[x,2] = c2
#     println("Done: ",ncodons[x])
# end

fac = 0.1/mutrate
gens = 10^(-2 + Int(floor(log10(fac))))
xmt = tmatX^gens;

xinit1 = orfvalsX[xProb,:].*Matrix(1.0*I,size(tmatX));
xinit1X = xinit1*xmt;

xinit2 = orfvalsX[xProb,:]'.*LowerTriangular(ones(size(tmatX)));

xinit3 = UpperTriangular(ones(size(tmatX))).*orfvalsX[xProb,:]';
setindex!.(Ref(xinit2), 0.0, eachindex(ncodons), eachindex(ncodons));
setindex!.(Ref(xinit3), 0.0, eachindex(ncodons), eachindex(ncodons));
xinit2X = xinit2*xmt;
xinit3X = xinit3*xmt;


# xinit1Y = UpperTriangular(xinit1X)
# setindex!.(Ref(xinit1Y),0.0,eachindex(ncodons),eachindex(ncodons));
maxcontri = [ncodons[findfirst(x .==maximum(x))] for x in eachcol(xinit1Y)];

# ptitle = replace((@sprintf "%0.0e" gens),"e+0" => "\\times10^{") * "}"
# ptitle = replace(ptitle,"1\\times" => "")

p1 = plot(ncodons[11:281],log10.(diag(xinit1X))[11:281],
    xlabel = "ORF length (codons)",
    ylabel = "log10(Probability)",
    label = "Length constant",
    title = latexstring(ptitle)* " Generations"
);
plot!(p1,ncodons[11:281],log10.(diag(xinit2X))[11:281], linecolor = :red, label = "Length extended");
plot!(p1,ncodons[11:281],log10.(diag(xinit3X))[11:281], linecolor = :blue, label = "Length truncated");

plot!(p1, legend = :bottom, legendfont = font(11,"Helvetica"),
    size = (width = cm2pt(13), height = cm2pt(11))
    )

savefig(p1,figdir*organism*"_Prob_"* (@sprintf "%0.0e" gens)*"gens.pdf")
