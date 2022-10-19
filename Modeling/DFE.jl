cd(Base.source_dir())
include("nucleotidefuncts.jl")

aasubdata = readdlm("BLOSUM62N.csv", ',', Int, '\n'; header = true)
aasubmat = zeros(21,21)
aasubmat[1:20,1:20] = aasubdata[1][1:20,1:20]
aasubmat[21,1:20] .= -100
aasubmat[1:20,21] .= -90
aanames = Dict(vcat([(aasubdata[2][i],i) for i in 1:20],("O",21)))
allcodons = kmers(3)

aasubmat[1:20,1:20] = (aasubmat[1:20,1:20] + aasubmat[1:20,1:20]')./2

aas = unique(gencode[:,3])
aaprob = [sum(nucprob.(gencode[gencode[:,3] .== x,1],gccontent)) for x in aas]

function nostop(codonset)
        codonset[codonset .∉ Ref(["TAA","TAG","TGA"])]
end

function dfecdn(codon1, codon2)
        mp = nucprob(codon1,gccontent)*mprob(codon1,codon2)
        fe = aasubmat[aanames[transnucs(codon1)],aanames[transnucs(codon2)]]
        return [-fe mp]
end

function dfedistgen(codonset1,codonset2)
        return vcat(reduce(vcat,[[dfecdn(y,x) for x in codonset2] for y in codonset1])...)
end

function combinedfe(dfeset)
        return reduce(vcat,[[i sum(dfeset[dfeset[:,1] .== i,2])] for i in unique(dfeset[:,1])])
end

function expdfe(dfeset)
        return sum(dfeset[:,1].*dfeset[:,2])/sum(dfeset[:,2])
end

function dfeframeX(cdn,frm)
        fset0 = frameXcodons(cdn,frm)
        # cset1 = singlemutant(cdn)
        cset1 = nostopcodons[nostopcodons .!= cdn]
        fset1 = hcat([frameXcodons(x,frm) for x in cset1]...)
        dfe1 = combinedfe(vcat(vcat([[dfecdn(fset0[y],x) for x in fset1[y,:]] for y in 1:20]...)...))
        return dfe1
end

nostopcodons = nostop(allcodons)

f1cdn = permutedims(hcat(frameXcodons.(nostopcodons,1)...))
f2cdn = permutedims(hcat(frameXcodons.(nostopcodons,2)...))

r1cdn = reverse.(comp.(f1cdn))
r2cdn = reverse.(comp.(f2cdn))

gccontent = 0.42

dfedist0 = combinedfe(dfedistgen(nostopcodons,allcodons))

frames = [-2,-1,1,2]

sdiff2 = (s1,s2) -> s1[s1 .∉ Ref(s2)]

dfeXframesehift = hcat([[dfeframeX(x,y) for x in nostopcodons] for y in frames]...)

[expdfe(vcat(dfeXframesehift[:,x]...)) for x in 1:4]