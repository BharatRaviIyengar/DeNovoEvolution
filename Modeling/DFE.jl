using Plots
using Measures
cd(Base.source_dir())
include("nucleotidefuncts.jl")

organism = "scer"

normalize = x -> x/sum(x)

cm2pt = (cm) -> 28.3465*cm
figdir = joinpath(Base.source_dir(),"../Manuscripts/Figures/M2_main/pdf/");
colors = ["#FFCC00","#5599FF","#D40000","#754473","#000000"];
lstyles = [:solid,:dash,:dot]

default(linecolor = :black, linewidth = 2, tickfont = font(10,"DejaVu Sans"), 
guidefont = font(13,"DejaVu Sans"),framestyle = :box, legend = false);

stopcodons = ["TAA","TAG","TGA"]
function nostop(codonset)
	codonset[codonset .∉ Ref(stopcodons)]
end

allcodons = kmers(3)
nostopcodons = nostop(allcodons);

nsub, nucsmbt = readdlm(organism*"_mutbias.txt",'\t',header = true);
if organism == "scer"
    mutrate = 1.7e-10
end

codonfreq = readdlm(joinpath(Base.source_dir(),organism*"_orf_codonfreq.txt"), '\t');

trimerfreq = readdlm(joinpath(Base.source_dir(),organism*"_intergenic_trimers.txt"), '\t');
trimerfreq[:,2] = normalize(trimerfreq[:,2])

PBMEC = 0 .+readdlm("PMBEC.txt")[2:21,2:21]
PBMEC2 = PBMEC .- 0.3;
setindex!.(Ref(PBMEC2), 0.0, 1:20, 1:20);
aanames = Dict(aas[i] => i for i in 1:20)
cdn2aanum = (cdn) -> aanames[transnucs(cdn)]

function mutdivcdn(gccontent)
    pmat,dmat = [zeros(61,60) for i = 1:2]
    for x in eachindex(nostopcodons)
        Z = nostopcodons[x]
		stZ = nostopcodons[nostopcodons .!= Z]
		pZ = nucprob(Z,gccontent)
        dmat[x,:] = [pZ*abs(PBMEC2[cdn2aanum(Z),cdn2aanum(y)]) for y in stZ]
        pmat[x,:] = [pZ*mprob(Z,y) for y in stZ]
    end
	dvec = vcat(dmat...)
	pvec = vcat(pmat...)/sum(pmat)
	dvx = unique(dvec)
    return [normalize(dvx) [sum(pvec[dvec .==q]) for q in dvx]]
end

function mutdivcdnX(t3)
	pmat,dmat = [zeros(61,60) for i = 1:2]
    for x in eachindex(nostopcodons)
        Z = nostopcodons[x]
		stZ = nostopcodons[nostopcodons .!= Z]
		pZ = nprob3(Z,t3)
        dmat[x,:] = [pZ*abs(PBMEC2[cdn2aanum(Z),cdn2aanum(y)]) for y in stZ]
        pmat[x,:] = [pZ*mprob(Z,y) for y in stZ]
    end
	dvec = vcat(dmat...)
	pvec = vcat(pmat...)/sum(pmat)
	dvx = unique(dvec)
    return [normalize(dvx) [sum(pvec[dvec .==q]) for q in dvx]]
end